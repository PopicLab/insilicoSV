from collections import defaultdict, Counter, namedtuple
from contextlib import closing
import logging
import os
import os.path
import re
from pysam import VariantFile
import copy
from intervaltree import IntervalTree, Interval

from insilicosv import utils
from insilicosv.utils import Region, Locus, if_not_none
from insilicosv.sv_defs import Operation, Transform, TransformType, BreakendRegion, VariantType, Syntax
from insilicosv.variant_set import get_vcf_header_infos

logger = logging.getLogger(__name__)

class StatsCollector:

    def __init__(self, chroms, chrom_lengths):
        assert len(chroms) == len(chrom_lengths)
        self.num_heterozygous = 0
        self.num_homozygous = 0
        self.placed_svs = 0
        self.total_svs = 0
        self.chroms = chroms
        self.chrom_lengths = chrom_lengths
        self.len_frags_chr = defaultdict(int)  # Lengths of altered fragments within chromosome
        self.sv_types = Counter()
        self.region_types = Counter()

    def get_info(self, svs):
        """
        collects all information for stats file after all edits are completed
        """
        self.total_svs = len(svs)
        
        self.min_region_len = None
        self.max_region_len = None

        tot_frags_len: int = 0
        tot_frags_count = 0

        for sv in svs:
            if not sv.is_placed():
                continue
            self.placed_svs += 1
            sv_type = sv.info.get('SVTYPE', 'UNKNOWN')
            self.sv_types[sv_type] += 1
            if (sv.roi is not None) and (sv.roi.kind != '_reference_'):
                self.region_types[sv.roi.kind] += 1
            assert sv.genotype is not None
            zygosity = sv.genotype[0] and sv.genotype[1] and (True if sv_type != VariantType.SNP else (
                    sv.replacement_seq[0] == sv.replacement_seq[1]))
            if zygosity:
                self.num_homozygous += 1
            else:
                self.num_heterozygous += 1

            for src_reg in set([operation.source_region for operation in sv.operations
                                if operation.source_region is not None]):
                self.len_frags_chr[src_reg.chrom] += src_reg.length()

                tot_frags_len += src_reg.length()
                tot_frags_count += 1

                if self.min_region_len is None or src_reg.length() < self.min_region_len:
                    self.min_region_len = src_reg.length()
                if self.max_region_len is None or src_reg.length() > self.max_region_len:
                    self.max_region_len = src_reg.length()
        self.avg_len = tot_frags_len // tot_frags_count if tot_frags_count != 0 else 0

    def export_data(self, fileout):
        """
        Exports all collected data to entered file
        fileout: Location to export stats file
        """
        def write_item(fout, name, item, prefix=""):
            fout.write("{}{}: {}\n".format(prefix, str(name), str(item)))

        os.makedirs(os.path.dirname(fileout), exist_ok=True)
        with open(fileout, "w") as fout:
            fout.write("===== Overview =====\n")
            write_item(fout, "SVs successfully simulated", str(self.placed_svs) + "/" + str(self.total_svs))
            for sv_type in self.sv_types:
                write_item(fout, sv_type, self.sv_types[sv_type], prefix="\t- ")
            write_item(fout, "Homozygous SVs", self.num_homozygous)
            write_item(fout, "Heterozygous SVs", self.num_heterozygous)
            write_item(fout, "Average length of impacted reference regions", self.avg_len)
            write_item(fout, "Min length of impacted reference region", self.min_region_len)
            write_item(fout, "Max length of impacted reference region", self.max_region_len)
            if self.region_types:
                write_item(fout, "Number of SVs per ROI types", self.region_types)
            for chrom, chrom_length in zip(self.chroms, self.chrom_lengths):
                fout.write(f"\n===== {chrom} =====\n")
                write_item(fout, "Length of sequence", chrom_length)
                write_item(fout, "Total impacted length of reference chromosome",
                           self.len_frags_chr[chrom])
# end: class StatsCollector

PafRecord = namedtuple('PafRecord', [
    "query_name",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target_name",
    "target_length",
    "target_start",
    "target_end",
    "residue_matches",
    "alignment_block_length",
    "mapping_quality",
    "tags",
])

class OutputWriter:

    def __init__(self, svs, recurrent_svs, reference, chrom_lengths, output_path, config):
        self.svs = svs
        self.svs = sorted(self.svs, key=operator.attrgetter('time_point'))
        self.recurrent_svs = recurrent_svs
        self.reference = reference
        self.chrom_lengths = chrom_lengths
        self.output_path = output_path
        self.config = config

    def output_haps(self):
        if self.config.get('output_no_haps', False):
            logger.warning('Skipping haps output')
            return
        for hap_index, hap_fa in enumerate(['sim.hapA.fa', 'sim.hapB.fa']):
            self.output_hap(os.path.join(self.output_path, hap_fa), hap_index, self.config.get('homozygous_only', False))

    def output_hap(self, hap_fa, hap_index, homozygous):
        paf_records = []
        hap_chrom_lengths = {}
        with open(hap_fa, 'w') as sim_fa:
            chrom2operations, target_regions, source_regions = self.get_chrom2operations(hap_index)
            for chrom, chrom_length in zip(self.reference.references, self.reference.lengths):
                hap_str = ['hapA', 'hapB'][hap_index]
                hap_chrom = f'{chrom}_{hap_str}'
                sim_fa.write(f'>{hap_chrom}\n')
                hap_chrom_pos = 0
                for overlapping_operations, sv_target, sv_source in zip(chrom2operations[chrom], target_regions[chrom], sv_sources[chrom]):
                    orig_seq = self.reference.fetch(
                                            reference=chrom,
                                            start=sv_source.start,
                                            end=sv_source.end)
                    # Tree used to keep track of the operations and different paf records
                    operation_tree = IntervalTree()

                    for operation in overlapping_operations:
                        # paf format: https://cran.r-project.org/web/packages/pafr/vignettes/Introduction_to_pafr.html
                        paf_mapq = 60
                        n_copies = operation.transform.n_copies
                        if operation.motif is not None:
                            # In the trEXP case the number of copies is encoded in the motif
                            n_copies = 1

                        for _ in range(n_copies):
                            if operation.transform_type == TransformType.DEL:
                                operation_start = operation.source_region.start
                                operation_end = operation.source_region.end
                                if 'side' in operation.op_info:
                                    # The recurrent operation was crossing over the target insertion point
                                    operation_length = operation_end - operation_start
                                    if operation.op_info['side'] == 0:
                                        operation_start = sv_source.start
                                        operation_end = operation_start + operation_length
                                    else:
                                        operation_end = sv_source.end
                                        operation_start = operation_end - operation_length
                                overlaps = operation_tree.overlap(operation_start, operation_end)
                                for overlap in overlaps:
                                    # This operation is included in a DEL so deleted
                                    operation_tree.remove(overlap)
                                    if overlap.begin < operation_start or overlap.end > operation_end:
                                        # We trim the interval and the associated PAF record to add them back to the tree
                                        if operation_start > overlap.begin:
                                            paf_trim = copy(overlap.data)
                                            if paf_trim.query_start == overlap.begin:
                                                paf_trim.query_start = overlap.begin
                                                paf_trim.query_length = operation_start - overlap.begin
                                            else:
                                                paf_trim.target_start = overlap.begin
                                                paf_trim.target_length = operation_start - overlap.begin
                                            operation_tree.add(Interval(begin=overlap.begin, end=operation_start, data=paf_trim))

                                        if operation_end < overlap.end:
                                            paf_trim = copy(overlap.data)
                                            if paf_trim.query_start == overlap.begin:
                                                paf_trim.query_start = operation_end
                                                paf_trim.query_length = overlap.end - operation_end
                                            else:
                                                paf_trim.target_start = operation_end
                                                paf_trim.target_end = overlap.end
                                            operation_tree.add(Interval(begin=operation_end, end=overlap.end, data=paf_trim))

                            if operation.transform_type == TransformType.DEL:
                                target_end = operation.source_region.end
                                length = operation.source_region.length()
                                if operation.motif is not None:
                                    length = len(operation.motif)
                                    target_end = operation.source_region.start + length
                                    # the rest of the tandem repeat region is kept
                                    seq = self.reference.fetch(
                                            reference=operation.source_region.chrom,
                                            start=target_end,
                                            end=operation.source_region.end)
                                paf_rec = PafRecord(
                                    query_name=hap_chrom, query_length=None,
                                    query_start=hap_chrom_pos, query_end=hap_chrom_pos,
                                    strand='+',
                                    target_name=operation.source_region.chrom,
                                    target_length=self.chrom_lengths[operation.source_region.chrom],
                                    target_start=operation.source_region.start,
                                    target_end=target_end,
                                    residue_matches=0, alignment_block_length=length,
                                    mapping_quality=paf_mapq,
                                    tags=f'cg:Z:{operation.source_region.length()}D')
                            else:
                                if operation.source_region is not None:
                                    if operation.motif is None:
                                        seq = self.reference.fetch(
                                            reference=operation.source_region.chrom,
                                            start=operation.source_region.start,
                                            end=operation.source_region.end)
                                    else:
                                        seq = operation.motif
                                    paf_rec = PafRecord(
                                        query_name=hap_chrom, query_length=None,
                                        query_start=hap_chrom_pos, query_end=hap_chrom_pos + len(seq),
                                        strand='+',
                                        target_name=operation.source_region.chrom,
                                        target_length=self.chrom_lengths[operation.source_region.chrom],
                                        target_start=operation.source_region.start,
                                        target_end=operation.source_region.end,
                                        residue_matches=len(seq), alignment_block_length=len(seq),
                                        mapping_quality=paf_mapq,
                                        tags=f'cg:Z:{len(seq)}M')
                                else:
                                    assert operation.novel_insertion_seq is not None
                                    seq = operation.novel_insertion_seq
                                    paf_rec = PafRecord(
                                        query_name=hap_chrom, query_length=None,
                                        query_start=hap_chrom_pos, query_end=hap_chrom_pos + len(seq),
                                        strand='+',
                                        target_name=operation.target_region.chrom,
                                        target_length=self.chrom_lengths[operation.target_region.chrom],
                                        target_start=operation.target_region.start,
                                        target_end=operation.target_region.end,
                                        residue_matches=0, alignment_block_length=len(seq),
                                        mapping_quality=paf_mapq,
                                        tags=f'cg:Z:{len(seq)}I')

                                if operation.transform_type == TransformType.INV:
                                    seq = utils.reverse_complement(seq)
                                    paf_rec = paf_rec._replace(strand='-')
                                if operation.transform.divergence_prob > 0:
                                    orig_seq = seq
                                    if operation.transform.replacement_seq is None or operation.transform.replacement_seq[hap_index] is None:
                                        haplotypes = operation.transform.replacement_seq
                                        if haplotypes is None:
                                            haplotypes = [None, None]
                                        if homozygous and hap_index == 1:
                                            replacement_seq = haplotypes[0]
                                        else:
                                            replacement_seq = utils.divergence(seq, operation.transform.divergence_prob)
                                        haplotypes[hap_index] = replacement_seq
                                        operation.transform = Transform(operation.transform_type,
                                                                        is_in_place=operation.transform.is_in_place,
                                                                        n_copies=operation.transform.n_copies,
                                                                        divergence_prob=operation.transform.divergence_prob,
                                                                        replacement_seq=haplotypes, orig_seq=orig_seq)
                                    seq = operation.transform.replacement_seq[hap_index]
                                    assert len(seq) == len(orig_seq)
                            if seq is not None:
                                sim_fa.write(seq)
                                hap_chrom_pos += len(seq)

                            # end: if operation.transform_type == TransformType.DEL:

                            assert paf_rec is not None

                            if operation.op_info and operation.op_info.get('SVID'):
                                sv_id = operation.op_info['SVID']
                                paf_rec = paf_rec._replace(tags=paf_rec.tags + f'\tsv:Z:{sv_id}')
                            paf_records.append(paf_rec)
                        # end: for _ in range(operation.transform.n_copies):
                    # end: for operation in overlapping_operations
                # end: for overlapping_operations in chrom2operations[chrom]

                hap_chrom_lengths[hap_chrom] = hap_chrom_pos
                sim_fa.write('\n')
            # end: for chrom, chrom_length in zip(...)
        # end: with open(hap_fa, 'w') as sim_fa
        if self.config.get('output_paf', False):
            hap_paf = hap_fa.replace('.fa', '.paf')

            new_paf_records = []
            for paf_rec_num, paf_rec in enumerate(paf_records):
                if (paf_rec.target_start == paf_rec.target_end and
                    paf_rec.query_start == paf_rec.query_end):
                    continue
                new_paf_records.append(paf_rec)

            with open(hap_paf, 'w') as sim_paf:
                for paf_rec in new_paf_records:
                    paf_rec = paf_rec._replace(query_length=hap_chrom_lengths[paf_rec.query_name])
                    sim_paf.write('\t'.join(map(str, paf_rec)) + '\n')

    # end: def output_hap(self, hap_fa, hap_index)

    def output_novel_adjacencies(self):
        if not self.config.get('output_adjacencies', False):
            logger.warning('Skipping adjacencies output')
            return
        novel_adjacencies = []
        for sv in self.recurrent_svs + self.svs:
            # Prevent duplicate adjacencies
            if sv.info['OP_TYPE'] in ['SNP', 'trEXP', 'trCON']: continue
            lhs, rhs = sv.info['GRAMMAR'].split('->')
            # Add symbols to represent the bases before and after the SV.
            lhs = ['PR'] + [symbol for symbol in lhs.strip() if
                              symbol not in [Syntax.ANCHOR_START, Syntax.ANCHOR_END, Syntax.DIVERGENCE]] + ['SU']
            rhs = ['PR'] + [symbol for symbol in rhs.strip() if
                              symbol not in [Syntax.ANCHOR_START, Syntax.ANCHOR_END, Syntax.DIVERGENCE]] + ['SU']

            # Get the original adjacencies and symbols
            breakends = {'PR': ['NA', sv.placement[0]], 'SU': [sv.placement[-1], 'NA']}
            lhs_adjacencies = []
            num_dispersion = 0
            prev_symbol = 'PR'
            for idx, symbol in enumerate(lhs[1:]):
                if symbol == Syntax.DISPERSION:
                    # to distinguish the different dispersions
                    symbol += str(num_dispersion)
                    num_dispersion += 1
                if symbol != 'SU':
                    breakends[symbol] = [sv.placement[idx], sv.placement[idx + 1]]
                # Adjacency between the end of the last symbol and the beginning of the current one.
                lhs_adjacencies.append([prev_symbol + '^t', symbol + '^h'])
                prev_symbol = symbol

            # Get the novel adjacencies after the SV placement
            prev_symbol = 'PR'
            num_dispersion = 0
            for symbol in rhs[1:]:
                if symbol == Syntax.DISPERSION:
                    # to distinguish the different dispersions
                    symbol += str(num_dispersion)
                    num_dispersion += 1
                elif symbol == Syntax.MULTIPLE_COPIES:
                    symbol = prev_symbol
                prev_orientation = '^t'
                prev_position = 1
                curr_position = 0
                curr_orientation = '^h'
                curr_symbol = symbol
                if prev_symbol.islower():
                    # Inversion
                    prev_orientation = '^h'
                    prev_symbol = prev_symbol.upper()
                    prev_position = 0
                if curr_symbol.islower():
                    curr_orientation = '^t'
                    curr_symbol = curr_symbol.upper()
                    curr_position = 1
                adjacency = [prev_symbol + prev_orientation, curr_symbol + curr_orientation]
                if adjacency not in lhs_adjacencies:
                    # Novel adjacency
                    locus_start = 'INS'
                    locus_end = 'INS'
                    if prev_symbol in breakends:
                        # Not a novel insertion
                        locus_start = breakends[prev_symbol][prev_position]
                    if curr_symbol in breakends:
                        # Not a novel insertion
                        locus_end = breakends[curr_symbol][curr_position]
                    genotype = str(int(sv.genotype[0])) + '|' + str(int(sv.genotype[1]))
                    record = [sv.info["GRAMMAR"], sv.sv_id, adjacency, locus_start, locus_end, genotype]
                    if record not in novel_adjacencies:
                        novel_adjacencies.append(record)
                prev_symbol = symbol
        adjacency_path = os.path.join(self.output_path, 'adjacencies.bed')
        with open(adjacency_path, 'w') as adjacency_file:
            for sv_grammar, sv_id, grammar, left_locus, right_locus, genotype in novel_adjacencies:
                chrom_start = pos_start = chrom_end = pos_end = pos_next_end = pos_next_start = 'INS'
                if left_locus != 'INS':
                    chrom_start = left_locus.chrom
                    pos_start = left_locus.pos
                    if '^t' in grammar[0]:
                        pos_start -= 1
                    pos_next_start = pos_start + 1
                if right_locus != 'INS':
                    chrom_end = right_locus.chrom
                    pos_end = right_locus.pos
                    if '^t' in grammar[1]:
                        pos_end -= 1
                    pos_next_end = pos_end + 1
                # Novel insertions are disregarded as they are not in the reference.
                if left_locus == 'INS' or right_locus == 'INS': continue
                record = [chrom_start, pos_start, pos_next_start, chrom_end, pos_end, pos_next_end,
                             '/'.join(grammar), sv_grammar, genotype, sv_id]
                adjacency_file.write('\t'.join(map(str, record)) + '\n')

    def get_chrom2operations(self, hap_index):
        """For each chromosome, make a list of operations targeting that chromosome,
        sorted along the chromosome.  Add identity operations for regions not touched
        by SVs."""
        merged_lookup = dict()
        chrom2operations = defaultdict(list)
        target_region2transform = dict(IntervalTree)
        for sv in recurrent_svs + self.svs:
            # start with recurrent SVs, regular SVs are ordered by time point
            assert sv.genotype is not None
            if sv.genotype[hap_index]:
                operations = []
                for operation in sv.operations:
                    operation.op_info['SVID'] = sv.sv_id
                    operation.op_info['recurrent'] = not sv.recurrent_freq is None
                    operation.op_info['time_point'] = sv.time_point
                    operations.append(operation)

                while operations:
                    operation = operations.pop()
                    assert operation.target_region is not None
                    chrom = operation.target_region.chrom
                    target_interval = Interval(operation.target_region.start, operation.end, data=operation)
                    overlap = target_region2transform[chrom].overlap(target_interval.begin, target_interval.end)
                    if operation.is_in_place:
                        if not operation.op_info['recurrent'] and any(not interval.data.op_info['recurrent'] for interval  in overlap):
                            # Check if there is already an inplace operation there, in which case they have to be the same
                            assert (operation.transform == target_region2transform[operation.target_region])
                            print(operation.transform, target_region2transform[operation.target_region])
                            continue

                    if operation.op_info is None:
                        operation.op_info = dict()
                    lookup_idx = len(chrom2operations[chrom])
                    merged_lookup[operation] = lookup_idx
                    chrom2operations[chrom].append([operation])

                    # Overlapping operations are stored in the same list to be applied together
                    for interval in overlap:
                        overlap_operation = interval.data
                        if not overlap_operation.op_info['recurrent'] or not overlap_operation in merged_lookup: continue
                        inter_start = max(interval.begin, target_interval.begin)
                        inter_end = min(interval.end, target_interval.end)

                        # Cut the recurrent operations to fit the target region
                        side = None
                        if target_interval.begin == target_interval.end:
                            # The recurrent operation is containing an insertion target, we decide if the right or left side will be included in the target.
                            side = np.argmin([target_interval.start - inter_start, inter_end - target_interval.end])

                        if not operation.op_info['recurrent'] and ((inter_start < interval.begin and inter_end == interval.end) or (side and side == 0)):
                            # the recurrent operation is on the right side, the operation outside the target region is added again to the operations to merge
                            outside_operation = utils.update_operation_groups(operation, interval.begin, target_interval.end,
                                                               target_interval.end, interval.end, merged_lookup,
                                                               target_region2transform, lookup_idx)
                            if side:
                                # The target is an insertion point and the recurrent operation is overlapping it
                                outside_operation.op_info['source_region'] = (operation_outside.source_region.start, operation_outside.source_region.end)
                                outside_operation.op_info['side'] = 1
                                operation.source_region.start = operation.source_region.end = target_interval.end
                            else:
                                operations.append(outside_operation)

                        elif not operation.op_info['recurrent'] and ((inter_start == interval.begin and inter_end < interval.end) or (side and side == 1)):
                            # the recurrent operation is on the left side, the operation outside the target region is added again to the operations to merge
                            outside_operation = utils.update_operation_groups(operation, target_interval.begin, interval.end,
                                                    interval.begin, target_interval.begin, merged_lookup,
                                                    target_region2transform, lookup_idx)
                            if side:
                                # The target is an insertion point and the recurrent operation is overlapping it
                                outside_operation.op_info['source_region'] = (operation_outside.source_region.start, operation_outside.source_region.end)
                                outside_operation.op_info['side'] = 0
                                operation.source_region.start = operation.source_region.end = target_interval.end
                            else:
                                operations.append(outside_operation)

                        else:
                            # The recurrent operation is included in the target
                            chrom2operations[chrom][lookup_idx].extend(chrom2operations[chrom][merged_lookup[overlap_operation]])
                            merged_lookup[overlap_operation] = lookup_idx
                            chrom2operations[chrom][merged_lookup[overlap_operation]] = []

                    target_region2transform[chrom].add(target_interval)
                    if not operation.is_in_place:
                        # The SV is not in place, reference the source region to keep track of past and future modifications
                        source_interval = Interval(operation.source_region.start, operation.source_region.end)
                        target_region2transform[chrom].overlap(source_interval.begin, source_interval.end)
                        for operation_overlap in overlap_source:
                            if not operation_overlap.op_info['recurrent']: continue
                            # All recurrent operations overlapping with the source and happening before have to impact the target
                            # If an operation happened on the source afterwards and is fully contained in the source, 50% chance to be on the target instead.
                            if operation.op_info['time_point']> operation_overlap.op_info['time_point']:
                                inter_start = max(operation_overlap.begin, source_interval.begin)
                                inter_end = min(operation_overlap.end, source_interval.end)
                                past_operation = Operation(transform=operation_overlap.transform,
                                                           source_breakend_region=BreakendRegion(0, 1),
                                                           placement=[Locus(chrom, inter_start),
                                                                      Locus(chrom, inter_end)])
                                chrom2operations[chrom][lookup_idx].append(past_operation)
                            elif (operation_overlap in merged_lookup) and (source_interval.begin < operation_overlap.begin < operation_overlap.end < source_interval.end):
                                # Ulterior operation unclaimed by an inplace operation and fully contained in the source
                                if random.randint(0, 1):
                                    chrom2operations[chrom][lookup_idx].append(
                                        chrom2operations[chrom][merged_lookup[operation_overlap]])
                                    # The overlap_operation is fully contained so keeping track of operation overlaps is enough.
                                    del merged_lookup[operation_overlap]


        # Clean up and order the operations by time point.
        for chrom in chrom2operations:
            chrom2operations[chrom] = [overlapping_operations for overlapping_operations in chrom2operations[chrom]
                                       if overlapping_operations]

            chrom2operations[chrom] = [sorted(operations, key=operator.attrgetter('time_point')) for operations in chrom2operations[chrom]]

        # Find the disjoint regions affected  by at least one operation.
        transformed_regions, source_regions = utils.get_transformed_regions(chrom2operations)
        for chrom, chrom_length in self.chrom_lengths.items():
            target_regions = sorted([Region(chrom=chrom, start=0, end=0)] + 
                                    [region for region in transformed_regions[chrome]] +
                                    [Region(chrom=chrom, start=chrom_length, end=chrom_length)])

            # Add identity operation between transformed regions to reconstruct the haplotype
            unchanged_regions = [Region(start=target_region1.end, end=target_region2.start, chrom=chrom) for target_region1, target_region2 in utils.pairwise(target_regions)
                if target_region2.start > target_region1.end]
            intersv_ops = [Operation(
                transform=Transform(transform_type=TransformType.IDENTITY,
                                    is_in_place=True),
                source_breakend_region=BreakendRegion(0,1),
                placement=[Locus(chrom, target_region.start),
                           Locus(chrom, target_region.end)])
                for target_region in unchanged_regions]

            transformed_regions[chrom].extend(unchanged_regions)
            chrom2operations[chrom].extend(intersv_ops)
            regions_order_idx = np.argsort(transformed_regions[chrom])
            transformed_regions[chrom] = sorted(transformed_regions[chrom])
            chrom2operations[chrom] = [chrom2operations[chrom][idx] for idx in regions_order_idx]
        return chrom2operations, transformed_regions, source_regions

    def output_novel_insertions(self):
        with open(os.path.join(self.output_path, 'sim.novel_insertions.fa'), 'w') as (
                sim_novel_insertions_fa):
            for sv in self.svs:
                for op_num, operation in enumerate(sv.operations):
                    if operation.novel_insertion_seq:
                        seq_id = (
                            f'{sv.sv_id}_{op_num}_{operation.target_region.chrom}_'
                            f'{operation.target_region.start}')
                        sim_novel_insertions_fa.write(f'>{seq_id}\n')
                        sim_novel_insertions_fa.write(f'{operation.novel_insertion_seq}\n')

    def output_vcf(self):
        vcf_path = os.path.join(self.output_path, 'sim.vcf')
        with open(vcf_path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrm, chrm_len in zip(self.reference.references,
                                      self.reference.lengths):
                vcf.write("##contig=<ID=%s,length=%d>\n" % (chrm, chrm_len))
            vcf.write("#%s\n" % "\t".join(["CHROM", "POS", "ID",
                                           "REF", "ALT", "QUAL", "FILTER", "INFO",
                                           "FORMAT", "SAMPLE"]))

        with closing(VariantFile(vcf_path)) as vcf_file:

            for vcf_header_info in get_vcf_header_infos():
                vcf_file.header.info.add(**vcf_header_info)

            vcf_file.header.formats.add('GT', number=1, type='String', description="Genotype")

            with closing(VariantFile(vcf_path, 'w', header=vcf_file.header)) as vcf_out_file:

                vcf_records: list[dict] = []

                for sv in self.svs:
                    vcf_records.extend(sv.to_vcf_records(self.config))
                assert not utils.has_duplicates(vcf_rec['id'] for vcf_rec in vcf_records)
                for vcf_rec in sorted(vcf_records, key=lambda rec: (rec['contig'], rec['start'])):
                    self.check_vcf_record(vcf_rec)
                    vcf_record = vcf_out_file.header.new_record(**vcf_rec)  # type: ignore
                    assert len(vcf_record.samples) == 1
                    vcf_record.samples[0].phased = True
                    vcf_out_file.write(vcf_record)
            # end: with closing(VariantFile(vcf_path, 'w', header=vcf_file.header)) as vcf_out_file
        # end: with closing(VariantFile(vcf_path)) as vcf_file
    # end: def output_vcf(self)

    def check_vcf_record(self, vcf_rec):
        assert set(vcf_rec.keys()) >= {'id', 'contig', 'start', 'stop', 'qual',
                                       'filter', 'alleles', 'samples'}
        assert not any(c.isspace() or c == ';' for c in vcf_rec['id'])
        info = vcf_rec.get('info')
        assert isinstance(info, (type(None), dict))
        for info_key, info_val in vcf_rec.get('info', {}).items():
            if isinstance(info_val, (str, list)):
                for val in utils.as_list(info_val):
                    assert re.search(r'[\s;,=]', val) is None

    def output_svops_bed(self):
        if not self.config.get('output_svops_bed', False):
            return

        with open(os.path.join(self.output_path, 'sim.svops.bed'), 'w') as out:
            for sv in self.svs:
                for operation in sv.operations:
                    op_str = '_'.join([sv.sv_id, sv.info['SVTYPE'], operation.transform_type.value])
                    if operation.op_info and 'symbol' in operation.op_info:
                        op_str += ('_' + operation.op_info['symbol'].name)

                    out.write('\t'.join(map(str,
                                            [operation.target_region.chrom,
                                             operation.target_region.start,
                                             max(operation.target_region.end,
                                                 operation.target_region.start+1),
                                             op_str]))
                              + '\n')

                    if (operation.source_region is not None and
                        operation.source_region != operation.target_region):
                        out.write('\t'.join(map(str,
                                                [operation.source_region.chrom,
                                                 operation.source_region.start,
                                                 max(operation.source_region.end,
                                                     operation.source_region.start),
                                                 op_str + '_src']))
                                  + '\n')

    def output_stats(self):
        stats_collector = StatsCollector(chroms=list(self.reference.references),
                                         chrom_lengths=list(self.reference.lengths))
        stats_collector.get_info(self.svs)
        stats_collector.export_data(os.path.join(self.output_path, 'sim.stats.txt'))

# end: class OutputWriter
