from collections import defaultdict, Counter, namedtuple
from contextlib import closing
import logging
import os
import os.path
import re
from pysam import VariantFile
import copy
from functools import cmp_to_key

from insilicosv import utils
from insilicosv.utils import Region, Locus, if_not_none
from insilicosv.sv_defs import Operation, Transform, TransformType, BreakendRegion, VariantType, Syntax, Breakend
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
            zygosity = sv.genotype[0] and sv.genotype[1] and ((sv_type not in [VariantType.SNP]) or (
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

    def __init__(self, svs, overlap_sv_regions, reference, chrom_lengths, output_path, allow_hap_overlap, config):
        self.svs = svs
        self.reference = reference
        self.chrom_lengths = chrom_lengths
        self.aneuploidy_chrom_lengths = copy.deepcopy(self.chrom_lengths)
        self.output_path = output_path
        self.config = config
        self.overlap_sv_regions = overlap_sv_regions
        self.homozygous_only = config.get('homozygous_only', False)
        self.allow_hap_overlap = allow_hap_overlap

    def output_haps(self):
        if self.config.get('output_no_haps', False):
            logger.warning('Skipping haps output')
            return
        haploids = ['sim.hapA.fa', 'sim.hapB.fa']
        if self.config.get('haploid', False):
            haploids = ['sim.fa']
        for hap_index, hap_fa in enumerate(haploids):
            self.output_hap(os.path.join(self.output_path, hap_fa), hap_index)

    def output_hap(self, hap_fa, hap_index):
        # paf format: https://cran.r-project.org/web/packages/pafr/vignettes/Introduction_to_pafr.html
        paf_records = []
        hap_chrom_lengths = {}
        paf_mapq = 60
        with open(hap_fa, 'w') as sim_fa:
            chrom2operations: dict[str, list[Operation]] = self.get_chrom2operations(hap_index)
            for chrom, operations in chrom2operations.items():
                # Add chromosome copies for aneuploidy, the chromosome name in the operation is added to the aneuploidy_chrom_lengths
                if chrom not in self.aneuploidy_chrom_lengths:
                    for overlapping_op in operations:
                        self.aneuploidy_chrom_lengths[chrom] = overlapping_op[-1].source_region.length()
                        # Aneuploidy can only be applied to regions starting at bp 0
                        if overlapping_op[-1].target_region.start > 0: break
            for chrom, chrom_length in self.aneuploidy_chrom_lengths.items():
                hap_str = ['hapA', 'hapB'][hap_index]
                hap_chrom = f'{chrom}_{hap_str}'
                hap_chrom_pos = 0
                first_seq = True

                for overlapping_operations in chrom2operations[chrom]:
                    # The non-overlapping operation if any is at the end
                    sv_region = overlapping_operations[-1].target_region
                    if overlapping_operations[-1].source_region:
                        sv_region = overlapping_operations[-1].source_region

                    seq = self.reference.fetch(
                        reference=sv_region.chrom,
                        start=sv_region.start,
                        end=sv_region.end)
                    position_shifts = []
                    length_shift = 0
                    # To change the position referential and find the start index of the operation in seq
                    for op_idx, operation in enumerate(overlapping_operations):
                        relative_position = sv_region.start
                        n_copies = operation.transform.n_copies[hap_index]
                        if operation.motif is not None:
                            # In the trEXP case the number of copies is encoded in the motif
                            n_copies = 1
                        target_start = operation.target_region.start
                        operation_start = operation_end = target_start
                        if operation.source_region:
                            operation_start = operation.source_region.start
                            # Represents the end of the operation for position_shifts
                            operation_end = operation.source_region.end
                        # Represents the length of the sequence affected by the operation
                        operation_length = operation_end - operation_start
                        if op_idx == len(overlapping_operations) - 1:
                            # The original SV is encompassing all the overlapped INS, its length has to be adapted
                            operation_length += length_shift

                        total_shift = 0
                        position_idx = 0
                        # Get the current position of the SV adapted after INS and DEL
                        for idx, (shift_type, shift_start, shift_length) in enumerate(position_shifts):
                            if shift_start < operation_start:
                                # the current operation starts after, shift accordingly
                                current_shift = shift_length
                                if shift_type == 'DEL':
                                    current_shift = -shift_length
                                total_shift += current_shift
                            else:
                                # The remaining shift intervals are after the current operation
                                position_idx = idx
                                break

                        # If INS or DEL, we adapt the shift
                        if operation.novel_insertion_seq:
                            position_shifts.insert(position_idx,
                                                   ('INS', operation_start, len(operation.novel_insertion_seq)))
                            length_shift += len(operation.novel_insertion_seq)
                        elif operation.transform_type == TransformType.DEL:
                            position_shifts.insert(position_idx,
                                                   ('DEL', operation_start, operation_end - operation_start))

                        # Position of the operation inside the sequence
                        start_in_region = operation_start - relative_position

                        if operation.transform_type == TransformType.DEL:
                            if operation_length <= 50:
                                # INDEL we report the original sequence
                                orig_seq = self.reference.fetch(
                                    reference=operation.source_region.chrom,
                                    start=operation_start,
                                    end=operation_end)
                                operation.transform = operation.transform.replace(orig_seq=orig_seq)

                            # Apply the DEL
                            seq = seq[:start_in_region + total_shift] + seq[
                                                                        start_in_region + operation_length + total_shift:]

                            paf_rec = PafRecord(
                                query_name=hap_chrom, query_length=None,
                                query_start=hap_chrom_pos, query_end=hap_chrom_pos,
                                strand='+',
                                target_name=operation.source_region.chrom,
                                target_length=self.chrom_lengths[operation.source_region.chrom],
                                target_start=operation_start,
                                target_end=operation_end,
                                residue_matches=0, alignment_block_length=operation_length,
                                mapping_quality=paf_mapq,
                                tags=f'cg:Z:{operation.source_region.length()}D')
                        else:
                            # The operation is not a DEL
                            # Modify the sequence
                            start_in_region = start_in_region + total_shift
                            modified_seq = seq[start_in_region:start_in_region + operation_length]
                            if operation.novel_insertion_seq:
                                modified_seq = operation.novel_insertion_seq
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
                            else:

                                if operation.motif:
                                    # TR
                                    modified_seq = operation.motif

                                strand = '+'
                                if operation.transform_type == TransformType.INV:
                                    modified_seq = utils.reverse_complement(modified_seq)
                                    strand = '-'

                                paf_rec = PafRecord(
                                    query_name=hap_chrom, query_length=None,
                                    query_start=hap_chrom_pos, query_end=hap_chrom_pos + len(modified_seq),
                                    strand=strand,
                                    target_name=operation.source_region.chrom,
                                    target_length=self.chrom_lengths[operation.source_region.chrom],
                                    target_start=operation_start,
                                    target_end=operation_end,
                                    residue_matches=len(modified_seq), alignment_block_length=len(modified_seq),
                                    mapping_quality=paf_mapq,
                                    tags=f'cg:Z:{len(modified_seq)}M')

                                # For a DUP or an mCNV apply the correct number of copies
                                modified_seq = modified_seq * n_copies

                            if operation.transform.divergence_prob > 0 or operation.transform.replacement_seq:
                                # There is a divergence
                                if operation.transform.replacement_seq is None or operation.transform.replacement_seq[hap_index] is None:
                                    # Insure the haplotypes respect the genotype specified and store it for writing in the VCF
                                    replacement_seq = utils.divergence(modified_seq,
                                                                       operation.transform.divergence_prob)
                                    if not operation.transform.replacement_seq:
                                        haplotypes = [replacement_seq if hap == hap_index else None for hap in [0, 1]]
                                    elif not self.homozygous_only and operation.transform.divergence_prob == 1:
                                        # The SNP can be homozygous or heterozygous with two different alleles
                                        haplotypes = [operation.transform.replacement_seq[0], replacement_seq]
                                    else:
                                        # The SNP/DIVERGENCE is homozygous
                                        haplotypes = [operation.transform.replacement_seq[0],
                                                      operation.transform.replacement_seq[0]]

                                    # Retain the replacement_seq and orig_seq for applying to other copies and to write in the VCF output
                                    operation.transform = operation.transform.replace(replacement_seq=haplotypes,
                                                                                      orig_seq=modified_seq)

                                modified_seq = operation.transform.replacement_seq[hap_index]

                            if operation.op_info and operation.op_info.get('SVID'):
                                sv_id = operation.op_info['SVID']
                                paf_rec = paf_rec._replace(tags=paf_rec.tags + f'\tsv:Z:{sv_id}')

                            seq = seq[:start_in_region] + modified_seq + seq[start_in_region + operation_length:]
                        if paf_rec:
                            paf_records.append(paf_rec)
                    hap_chrom_pos += len(seq)
                    if seq:
                        if first_seq:
                            # Ensure we only write the chromosome_hap name if it is not empty (might cause issues when reading the file)
                            sim_fa.write(f'>{hap_chrom}\n')
                            first_seq = False
                        sim_fa.write(seq)
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
        for sv in self.svs:
            breakends = sorted(list(sv.placement.keys()))
            placement = [sv.placement[breakend] for breakend in breakends]
            # Prevent duplicate adjacencies
            if sv.info['OP_TYPE'] in ['SNP', 'trEXP', 'trCON']: continue
            lhs, rhs = sv.info['GRAMMAR'].split('->')
            # Add symbols to represent the bases before and after the SV.
            lhs = ['PR'] + [symbol for symbol in lhs.strip() if
                            symbol not in [Syntax.ANCHOR_START, Syntax.ANCHOR_END, Syntax.DIVERGENCE]] + ['SU']
            rhs = ['PR'] + [symbol for symbol in rhs.strip() if
                            symbol not in [Syntax.ANCHOR_START, Syntax.ANCHOR_END, Syntax.DIVERGENCE]] + ['SU']

            # Get the original adjacencies and symbols
            breakends = {'PR': ['NA', placement[0]], 'SU': [placement[-1], 'NA']}
            lhs_adjacencies = []
            num_dispersion = 0
            prev_symbol = 'PR'
            for idx, symbol in enumerate(lhs[1:]):
                if symbol == Syntax.DISPERSION:
                    # to distinguish the different dispersions
                    symbol += str(num_dispersion)
                    num_dispersion += 1
                if symbol != 'SU':
                    breakends[symbol] = [placement[idx], placement[idx + 1]]
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

    def group_overlap_operations(self, overlap_sv, is_placed, region_to_overlap_start,
                                 region_to_overlap_end, hap_index):
        overlapping_region = []
        placed_overlap_sv = []
        # Group overlapping operations together and potentially slice them if they are partially overlapping
        for overlap in overlap_sv:
            sv = overlap.data.sv

            if sv.genotype[hap_index] and not (overlap.data.end == region_to_overlap_start
                                               or region_to_overlap_end == overlap.data.start):
                # We do not consider SVs at the boundaries as overlapping
                overlapping_operation = sv.operations[0]
                if (overlap.data.start < region_to_overlap_start < overlap.data.end or
                        overlap.data.start < region_to_overlap_end < overlap.data.end):
                    # A non-contained overlapping DEL
                    overlapping_start = max(overlap.data.start, region_to_overlap_start)
                    overlapping_end = min(overlap.data.end, region_to_overlap_end)

                    # Add potential remaining left and right operations
                    if overlapping_start > overlap.data.start:
                        left_sv = copy.deepcopy(sv)
                        placement = {
                            Breakend(0): Locus(pos=overlap.data.start, chrom=overlapping_operation.target_region.chrom),
                            Breakend(1): Locus(pos=region_to_overlap_start,
                                               chrom=overlapping_operation.target_region.chrom)}

                        left_sv.operations[0].update_placement(placement)
                        self.overlap_sv_regions.add_region(Region(chrom=overlapping_operation.target_region.chrom,
                                                                  start=overlap.data.start,
                                                                  end=region_to_overlap_start),
                                                           sv=left_sv)

                    if overlapping_end < overlap.data.end:
                        right_sv = copy.deepcopy(sv)
                        placement = {
                            Breakend(0): Locus(pos=region_to_overlap_end,
                                               chrom=overlapping_operation.target_region.chrom),
                            Breakend(1): Locus(pos=overlap.data.end, chrom=overlapping_operation.target_region.chrom)}

                        right_sv.operations[0].update_placement(placement)

                        self.overlap_sv_regions.add_region(
                            Region(chrom=overlapping_operation.target_region.chrom,
                                   start=region_to_overlap_end,
                                   end=overlap.data.end),
                            sv=right_sv)

                    overlapping_operation.placement = {
                        Breakend(0): Locus(pos=overlapping_start, chrom=overlapping_operation.target_region.chrom),
                        Breakend(1): Locus(pos=overlapping_end, chrom=overlapping_operation.target_region.chrom)}

                    # Keep the information on the original position
                    overlapping_operation.orig_start = overlap.data.start
                    overlapping_operation.orig_end = overlap.data.end

                overlapping_region.append(overlapping_operation)
                if is_placed:
                    placed_overlap_sv.append(overlap)
        return overlapping_region, placed_overlap_sv

    def get_chrom2operations(self, hap_index):
        """For each chromosome, make a list of operations targeting that chromosome,
        sorted along the chromosome.  Add identity operations for regions not touched
        by SVs."""
        chrom2operations = defaultdict(list)
        target_region2transform = {}
        placed_overlap_sv = []
        # If not self.allow_hap_overlap, there is only one tree per chrom for the sv overlap else three
        hap_id_overlap = 0
        if self.allow_hap_overlap:
            hap_id_overlap = hap_index

        for sv in self.svs:
            # If the SV is on the other haplotype or overlapping it will be treated later
            if not sv.genotype[hap_index] or sv.allow_sv_overlap: continue

            for operation in sv.operations:
                assert operation.target_region is not None
                if operation.is_in_place:
                    if operation.target_region in target_region2transform:
                        assert (operation.transform ==
                                target_region2transform[operation.target_region])
                        continue
                    target_region2transform[operation.target_region] = operation.transform

                if operation.op_info is None:
                    operation.op_info = dict()
                operation.op_info['SVID'] = sv.sv_id
                if operation.transform.divergence_prob > 0:
                    operation.genotype = sv.genotype

                overlap_sv = self.overlap_sv_regions.chrom2itree[operation.target_region.chrom][hap_id_overlap].overlap(
                    operation.target_region.start - 0.1,
                    operation.target_region.end + 0.1)
                region_to_overlap_start = operation.target_region.start
                region_to_overlap_end = operation.target_region.end
                is_placed = True
                overlapping_region, placed_overlap_sv_cpt = self.group_overlap_operations(overlap_sv, is_placed,
                                                                                          region_to_overlap_start,
                                                                                          region_to_overlap_end,
                                                                                          hap_index)
                placed_overlap_sv += placed_overlap_sv_cpt

                if operation.source_region and operation.source_region != operation.target_region:
                    # In the case of a duplication, a change in the source has to be reflected on the target
                    overlap_sv = self.overlap_sv_regions.chrom2itree[operation.target_region.chrom][
                        hap_id_overlap].overlap(operation.source_region.start,
                                                operation.source_region.end)
                    region_to_overlap_start = operation.source_region.start
                    region_to_overlap_end = operation.source_region.end
                    # The overlapping SV still has to be placed on the source region
                    is_placed = False
                    overlapping_region_cpt, placed_overlap_sv_cpt = self.group_overlap_operations(overlap_sv, is_placed,
                                                                                                  region_to_overlap_start,
                                                                                                  region_to_overlap_end,
                                                                                                  hap_index)
                    placed_overlap_sv += placed_overlap_sv_cpt
                    overlapping_region += overlapping_region_cpt

                # Sort the overlapping operations to ensure determinism of the SNP sequences
                overlapping_region.sort(key=lambda op: op.target_region.start)
                overlapping_region.append(operation)
                chrom2operations[operation.target_region.chrom].append(overlapping_region)

        for chrom, trees in self.overlap_sv_regions.chrom2itree.items():
            tree = trees[hap_id_overlap]
            for overlap_regions in tree:
                sv = overlap_regions.data.sv
                if sv.genotype[hap_index] and not overlap_regions in placed_overlap_sv:
                    chrom2operations[chrom].append([sv.operations[0]])

        for chrom, chrom_length in self.chrom_lengths.items():
            target_regions = sorted([Region(chrom=chrom, start=0, end=0)] +
                                    [operations[-1].target_region for operations in chrom2operations[chrom]] +
                                    [Region(chrom=chrom, start=chrom_length, end=chrom_length)])

            intersv_ops = [[Operation(
                transform=Transform(transform_type=TransformType.IDENTITY,
                                    is_in_place=True),
                op_info={'SVID': 'sv' + str(idx_id + len(self.svs))},
                source_breakend_region=BreakendRegion(0,1),
                placement={Breakend(0): Locus(chrom, target_region1.end),
                           Breakend(1): Locus(chrom, target_region2.start)})]
                for idx_id, (target_region1, target_region2) in enumerate(utils.pairwise(target_regions))
                if target_region2.start > target_region1.end]

            chrom2operations[chrom].extend(intersv_ops)

            def compare_operations(operations1, operations2):
                # The last operation of the list is the non-overlapping one if any, if not they will have the same target
                operation1 = operations1[-1]
                operation2 = operations2[-1]

                # Primary sort: target_region_start
                if operation1.target_region.start != operation2.target_region.start:
                    return operation1.target_region.start - operation2.target_region.start

                # Secondary sort: target_region_end if at least one of them is not an insertion target
                if operation1.target_region.end != operation2.target_region.end:
                    return operation1.target_region.end - operation2.target_region.end

                # Tertiary sort: the operation comes from the same SV, the insertion order if provided has to be used
                if (operation1.target_insertion_order and operation2.target_insertion_order and
                        operation1.target_insertion_order[0] == operation2.target_insertion_order[0]):
                    return operation1.target_insertion_order[1] - operation2.target_insertion_order[1]

                # Quaternary case, different SVs, the insertion targets are inserted by proximity to the source if any (so DUPs are not separated from their copies)
                if operation1.source_region:
                    # Targets belong to the same chrom by design not source regions
                    if (operation1.source_region.chrom == operation1.target_region.chrom and
                            operation1.source_region.end == operation1.target_region.start):
                        # operation1 is inserted right next to its source, on the right, it has to be first
                        return -1
                    elif (operation1.source_region.chrom == operation1.target_region.chrom and
                          operation1.target_region.end == operation1.source_region.start):
                        # operation1 is inserted right next to its source, on the left, it has to be second
                        return 1

                if operation2.source_region:
                    if (operation2.source_region.chrom == operation2.target_region.chrom and
                            operation2.source_region.end == operation2.target_region.start):
                        # operation2 is inserted right next to its source, on the right, it has to be first
                        return 1
                    elif (operation2.source_region.chrom == operation2.target_region.chrom and
                          operation2.target_region.end == operation2.source_region.start):
                        # operation2 is inserted right next to its source, on the left, it has to be second
                        return -1
                # None of the operations are inserted at a position adjacent to the source, the order is arbitrary, we choose the svid
                return 1 if operation1.op_info['SVID'] > operation2.op_info['SVID'] else -1

            chrom2operations[chrom].sort(key=cmp_to_key(compare_operations))
        return chrom2operations

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

    def output_divergence(self):
        with open(os.path.join(self.output_path, 'sim.divergence.fa'), 'w') as (
                sim_novel_insertions_fa):
            for sv in self.svs:
                for op_num, operation in enumerate(sv.operations):
                    if 0 < operation.transform.divergence_prob < 1 and operation.transform.replacement_seq:
                        for hap_idx, hap in enumerate(['hapA', 'hapB']):
                            if not operation.transform.replacement_seq[hap_idx]: continue
                            seq_id = (
                                f'{sv.sv_id}_{op_num}_{operation.target_region.chrom}_'
                                f'{operation.target_region.start}_{hap}')
                            sim_novel_insertions_fa.write(f'>{seq_id}\n')
                            sim_novel_insertions_fa.write(f'{operation.transform.replacement_seq[hap_idx]}\n')

    def output_vcf(self):
        vcf_path = os.path.join(self.output_path, 'sim.vcf')
        with open(vcf_path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrm, chrm_len in self.aneuploidy_chrom_lengths.items():
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
                                                 operation.target_region.start + 1),
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