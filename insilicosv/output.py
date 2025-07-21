from collections import defaultdict, Counter, namedtuple
from contextlib import closing
import logging
import os
import os.path
import re
from pysam import VariantFile
from intervaltree import IntervalTree, Interval
import operator
import numpy as np
import random

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
        with open(hap_fa, 'w') as sim_fa:
            chrom2operations, source_regions = self.get_chrom2operations(hap_index)
            for chrom, chrom_length in zip(self.reference.references, self.reference.lengths):
                hap_str = ['hapA', 'hapB'][hap_index]
                hap_chrom = f'{chrom}_{hap_str}'
                sim_fa.write(f'>{hap_chrom}\n')
                for overlapping_operations, sv_source in zip(chrom2operations[chrom], source_regions[chrom]):
                    print(overlapping_operations)
                    seq = self.reference.fetch(
                                            reference=sv_source.chrom,
                                            start=sv_source.start,
                                            end=sv_source.end)
                    position_shifts = []
                    for operation in overlapping_operations:
                        n_copies = operation.transform.n_copies
                        if operation.motif is not None:
                            # In the trEXP case the number of copies is encoded in the motif
                            n_copies = 1
                        operation_start = operation_end = operation.target_region.start
                        if operation.source_region:
                            operation_start = operation.source_region.start
                            # Represents the end of the operation for position_shifts
                            operation_end = operation.source_region.end
                            if operation.source_region.overlap_position:
                                operation_start = operation_start + operation.source_region.overlap_position
                                operation_end = operation_start + operation.origin_length
                        # Represents the length of the sequence affected by the operation
                        operation_length = operation_end - operation_start + 1

                        total_shift = 0

                        if operation.transform_type == TransformType.DEL:
                            if 'side' in operation.op_info:
                                # The recurrent operation was crossing over the target insertion point
                                if operation.op_info['side'] == 0:
                                    operation_start = sv_source.start
                                    operation_end = operation_start + operation_length - 1
                                else:
                                    operation_end = sv_source.end
                                    operation_start = operation_end - operation_length + 1

                            # Get the current position of the SV adapted after INS and DEL
                            idx = 0
                            if position_shifts:
                                while position_shifts:
                                    shift_type, shift_start, shift_length = position_shifts[idx]
                                    if shift_start <= operation_start <= shift_start + shift_length - 1:
                                        # Overlapping DELs
                                        if operation_end <= shift_start + shift_length - 1:
                                            # The new DEL is included
                                            break
                                        operation_start = shift_start + shift_length - 1
                                        operation_length -= shift_length
                                    elif shift_start + shift_length - 1 < operation_start:
                                        # the current operation starts after, shift accordingly
                                        current_shift = shift_length
                                        if shift_type == 'DEL':
                                            current_shift = -shift_length
                                        total_shift += current_shift
                                    elif shift_start <= operation_end:
                                        if shift_type == 'DEL':
                                            # The rest has already been deleted, if INS we remove part of the INS sequence.
                                            operation_length -= min(shift_length, operation_end - shift_start + 1)
                                            shift_length -= operation_end - shift_start
                                        else:
                                            if operation_end > shift_start + shift_length:
                                                operation_end -= shift_length
                                                shift_length = 0
                                            else:
                                                shift_length -= operation_end - shift_start
                                                operation_start = shift_start
                                        if shift_start + shift_length > operation_end:
                                            # The shift operation has been truncated
                                            position_shifts.insert(idx, (shift_type, operation_end, shift_length))
                                        else:
                                            idx -= 1
                                    else:
                                        break
                                    idx += 1

                            # The remaining shift intervals are after the current operation
                            position_shifts.insert(idx, ('DEL', operation_start, operation_end - operation_start))
                            start_in_region = operation_start - sv_source.start
                            seq = seq[:start_in_region+total_shift] + seq[start_in_region+operation_length+total_shift:]
                        else:
                            # The operation is not a DEL
                            idx = 0
                            for shift_type, shift_start, shift_end, shift_length in position_shifts:
                                multiplier = -(shift_type == 'DEL')
                                idx += 1
                                if shift_start <= operation_start:
                                    if operation_end <= shift_start + shift_length - 1:
                                        total_shift += multiplier * (operation_end - shift_start + 1)
                                    else:
                                        total_shift += multiplier * shift_length
                                else:
                                    break
                            if operation.novel_insertion_seq:
                                # If insertion, it will cause a shift
                                position_shifts.insert(idx, ('INS', operation_start, len(operation.novel_insertion_seq)))

                            # Modify the sequence
                            operation_start = operation_start + total_shift - sv_source.start
                            modified_seq = seq[operation_start:operation_start+operation_length - 1]

                            if operation.novel_insertion_seq:
                                modified_seq = operation.novel_insertion_seq
                            if operation.transform_type == TransformType.INV:
                                modified_seq = utils.reverse_complement(modified_seq)
                            if n_copies > 1:
                                modified_seq = modified_seq * n_copies

                            if operation.transform.divergence_prob > 0:
                                # There is a SNP
                                if operation.transform.replacement_seq is None:
                                    # Insure the haplotypes respect the genotype specified and store it for writing in the VCF
                                    replacement_seq = utils.divergence(modified_seq, operation.transform.divergence_prob)
                                    if homozygous:
                                        haplotypes = [replacement_seq, replacement_seq]
                                    else:
                                        haplotypes = [modified_seq, modified_seq]
                                        haplotypes[hap_index] = replacement_seq
                                    operation.transform = Transform(operation.transform_type,
                                                                is_in_place=operation.transform.is_in_place,
                                                                n_copies=operation.transform.n_copies,
                                                                divergence_prob=operation.transform.divergence_prob,
                                                                replacement_seq=haplotypes, orig_seq=seq)
                                modified_seq = operation.transform.replacement_seq[hap_index]

                            seq = seq[:operation_start] + modified_seq + seq[operation_start + operation_length:]
                    print('SEQ', seq)
                    if seq:
                        sim_fa.write(seq)
                sim_fa.write('\n')


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

    @staticmethod
    def update_operation_groups(overlap_operation, overlap_start, overlap_end, target_operation, target_start, target_end, merged_lookup,
                                target_region2transform, chrom2operations, lookup_idx):
        # Split a recurrent operation into two operations, one inside and one outside the target interval, update the tree and lists of operations
        chrom = overlap_operation.target_region.chrom
        if target_start != target_end:
            # Inplace target operation
            inside_start = max(overlap_start, target_start)
            inside_end = min(overlap_end, target_end)
            inside_length = inside_end - inside_start

            operation_inside = Operation(transform=overlap_operation.transform,
                                         source_breakend_region=BreakendRegion(0, 1),
                                         placement=[Locus(chrom, inside_start),
                                                    Locus(chrom, inside_end)],
                                         op_id=overlap_operation.op_id + 'in',
                                         time_point=overlap_operation.time_point,
                                         recurrent=True)
        else:
            # Insertion target
            if overlap_operation.overlap_position:
                # The recurrent operation was placed inside the interval
                inside_position = overlap_operation.overlap_position
                inside_length = min(target_operation.origin_length - inside_position, overlap_operation.origin_length)
            else:
                # The recurrent SV starts before the operation and overlaps it.
                inside_position = 0
                inside_length = min(overlap_start + overlap_operation.origin_length - target_start, target_operation.origin_length)
            target_operation.origin_length -= inside_length
            operation_inside = Operation(transform=overlap_operation.transform,
                                         source_breakend_region=BreakendRegion(0, 1),
                                         placement=[Locus(chrom, target_start),
                                                    Locus(chrom, target_start)],
                                         overlap_position=inside_position,
                                         origin_length=inside_length,
                                         time_point=overlap_operation.time_point,
                                         op_id=overlap_operation.op_id + 'in',
                                         recurrent=True)

        operations_outside = []
        left_length = 0
        if overlap_start < target_start:
            # The overlap started on the left of the target
            left_length = target_start - overlap_start
            aux_outside = Operation(transform=overlap_operation.transform,
                                      source_breakend_region=BreakendRegion(0, 1),
                                      placement=[Locus(chrom, overlap_start),
                                                 Locus(chrom, target_start)],
                                      op_id=overlap_operation.op_id + 'out',
                                      time_point=overlap_operation.time_point,
                                     recurrent=True)
            operations_outside.append(aux_outside)
            # Update the tree to detect overlap with the outside operations
            target_region2transform[chrom].add(
                Interval(begin=overlap_start - 0.1, end=outside_end + 0.1, data=aux_outside))

        if overlap_end > target_end or inside_length + left_length < overlap_operation.origin_length:
            # The overlap started on the right or inside the target
            remaining_length = overlap_operation.origin_length - inside_length - left_length
            aux_outside = Operation(transform=overlap_operation.transform,
                                                source_breakend_region=BreakendRegion(0, 1),
                                                placement=[Locus(chrom, target_end + 1),
                                                           Locus(chrom, target_end + remaining_length)],
                                                op_id=overlap_operation.op_id + 'out',
                                                time_point=overlap_operation.time_point,
                                                recurrent=True)
            operations_outside.append(aux_outside)
            # Update the tree to detect overlap with the outside operations
            target_region2transform[chrom].add(
                Interval(begin=target_end + 1 - 0.1, end=target_end + remaining_length + 0.1, data=aux_outside))

        # Update the stored operations
        chrom2operations[chrom][merged_lookup[overlap_operation.op_id]].remove(overlap_operation)
        del merged_lookup[overlap_operation.op_id]
        # Merge the inside operation
        chrom2operations[chrom][lookup_idx].append(operation_inside)
        merged_lookup[operation_inside.op_id] = lookup_idx

        return operations_outside

    def get_chrom2operations(self, hap_index):
        """For each chromosome, make a list of operations targeting that chromosome,
        sorted along the chromosome.  Add identity operations for regions not touched
        by SVs."""
        merged_lookup = dict()
        chrom2operations = defaultdict(list)
        target_region2transform = defaultdict(IntervalTree)
        for sv in self.recurrent_svs + self.svs:
            # start with recurrent SVs, regular SVs are ordered by time point
            print('SV', sv, sv.genotype[hap_index])
            if sv.genotype[hap_index]:
                operations = []
                for op_id, operation in enumerate(sv.operations):
                    operation.op_info['SVID'] = sv.sv_id
                    operation.recurrent = sv.overlap_sv
                    operation.time_point = sv.time_point
                    operation.op_id = str(sv.sv_id) + '_' + str(op_id)
                    operations.append(operation)
                print('OPERATIONS', operations)
                while operations:
                    operation = operations.pop()
                    chrom = operation.target_region.chrom
                    target_start = operation.target_region.start
                    target_end = operation.target_region.end
                    # Padding of the interval to be able to include points as intervals.
                    target_interval = Interval(begin=target_start-0.1, end=target_end+0.1, data=operation)

                    # Initialize the length of the target available for overlap
                    if not operation.origin_length:
                        if operation.source_region:
                            length = operation.source_region.length()
                        else:
                            length = len(operation.novel_insertion_seq)
                        operation.origin_length = length

                    overlap = target_region2transform[chrom].overlap(target_interval.begin, target_interval.end)
                    if operation.is_in_place:
                        if target_region2transform[operation.target_region] and (not operation.recurrent and any(not interval.data.recurrent for interval  in overlap)):
                            # Check if there is already an inplace operation there, in which case they have to be the same
                            assert (operation.transform in target_region2transform[operation.target_region])
                            continue

                    if operation.op_info is None:
                        operation.op_info = dict()
                    lookup_idx = len(chrom2operations[chrom])
                    merged_lookup[operation.op_id] = lookup_idx
                    chrom2operations[chrom].append([operation])

                    # Overlapping operations are stored in the same list to be applied together
                    for overlap_interval in overlap:
                        overlap_operation = overlap_interval.data
                        overlap_start = overlap_interval.data.target_region.start
                        overlap_end = overlap_interval.data.target_region.end
                        # The overlap operation is not recurrent or has been merged with a non-recurrent SV
                        if not overlap_operation.recurrent or not overlap_operation.op_id in merged_lookup: continue

                        # Merge the part of the overlap_operation overlapping the current target, append to the operations the remaining part
                        outside_operations = self.update_operation_groups(overlap_operation, overlap_start, overlap_end,
                                                                          operation,
                                                                         target_start, target_end, merged_lookup,
                                                                         target_region2transform, chrom2operations,
                                                                         lookup_idx)
                        if outside_operations:
                            operations += outside_operations

                    target_region2transform[chrom].add(target_interval)
                    if not operation.is_in_place and operation.source_region:
                        # The SV is not in place and not a novel insertion, reference the source region to keep track of past and future modifications
                        source_interval = Interval(operation.source_region.start, operation.source_region.end)
                        overlap_source = target_region2transform[chrom].overlap(source_interval.begin, source_interval.end)
                        for interval_overlap in overlap_source:
                            operation_overlap = interval_overlap.data
                            if not operation_overlap.recurrent: continue
                            # All recurrent operations overlapping with the source and happening before have to impact the target
                            if operation.time_point> interval_overlap.time_point:
                                inside_start = max(interval_overlap.begin, source_interval.begin)
                                inside_end = min(interval_overlap.end, source_interval.end)
                                past_operation = Operation(transform=interval_overlap.transform,
                                                           source_breakend_region=BreakendRegion(0, 1),
                                                           placement=[Locus(chrom, inside_start),
                                                                      Locus(chrom, inside_end)],
                                                           time_point=interval_overlap.time_point)
                                chrom2operations[chrom][lookup_idx].append(past_operation)

        # Clean up and order the operations by time point.
        for chrom in chrom2operations:
            chrom2operations[chrom] = [overlapping_operations for overlapping_operations in chrom2operations[chrom]
                                       if overlapping_operations]
            chrom2operations[chrom] = [sorted(operations, key=lambda ope: ope.time_point) for operations in chrom2operations[chrom]]

        # Find the disjoint regions affected  by at least one operation.
        transformed_regions, source_regions = utils.get_transformed_regions(chrom2operations)
        for chrom, chrom_length in self.chrom_lengths.items():
            target_regions = sorted([Region(chrom=chrom, start=0, end=0)] + 
                                    [region for region in transformed_regions[chrom]] +
                                    [Region(chrom=chrom, start=chrom_length, end=chrom_length)])

            # Add identity operation between transformed regions to reconstruct the haplotype
            unchanged_regions = [Region(start=target_region1.end, end=target_region2.start, chrom=chrom) for target_region1, target_region2 in utils.pairwise(target_regions)
                if target_region2.start > target_region1.end]
            intersv_ops = [[Operation(
                transform=Transform(transform_type=TransformType.IDENTITY,
                                    is_in_place=True),
                source_breakend_region=BreakendRegion(0,1),
                placement=[Locus(chrom, target_region.start),
                           Locus(chrom, target_region.end)])]
                for target_region in unchanged_regions]

            transformed_regions[chrom].extend(unchanged_regions)
            source_regions[chrom].extend(unchanged_regions)
            chrom2operations[chrom].extend(intersv_ops)
            # The operations are directly written in the order of the target locations
            regions_order_idx = np.argsort(transformed_regions[chrom])
            chrom2operations[chrom] = [chrom2operations[chrom][idx] for idx in regions_order_idx]
            source_regions[chrom] = [source_regions[chrom][idx] for idx in regions_order_idx]
        return chrom2operations, source_regions

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
