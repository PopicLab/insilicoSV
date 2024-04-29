from collections import defaultdict
import argparse
import copy
import itertools
import os
import random

import numpy as np
from insilicosv.constants import Symbols
from insilicosv.processing import FormatterIO
from collections import defaultdict
from pysam import VariantFile


class NestedDict(defaultdict):
    def __call__(self):
        return NestedDict(self.default_factory)


def is_overlapping(event_ranges, addition, called_from_helper=False, strictly_partial=False):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored
    # 'strictly_partial' a toggle to return True if the interval overlaps a stored event, without being equal to it
    for event in event_ranges:
        if event[1] > addition[0] and event[0] < addition[1]:
            if called_from_helper:
                raise Exception("Overlap between {} and {}".format(event[0:2], addition[0:2]))
            else:
                return True if not strictly_partial else event != addition

    return False


def fail_if_any_overlapping(arr):
    # will raise Exception if any overlap between intervals is found
    # arr: list of tuples
    for x, ele in enumerate(arr):
        is_overlapping(arr[:x], ele, called_from_helper=True)


def validate_symbols(source, target):
    '''
    Ensures that source has unique symbols other than dispersion events - specifically used for source sequence
    Raises error if number of dispersions between source and target is not the same

    source: tuple, source sequence
    target: tuple, target sequence
    '''
    present = dict()
    for symbol in source:
        if symbol in present:
            raise Exception("Source transformation {} does not have unique symbols!".format(source))
        elif symbol != Symbols.DIS.value:   # exclude dispersion events because they will always appear the same for user inputs
            present[symbol] = True

        if Symbols.DUP.value in symbol:
            raise Exception("Duplication marking (') not allowed in source sequence {}".format(source))
        if any(c.islower() for c in symbol):
            raise Exception("Only original symbols may appear in source sequence, and they must also be uppercase: {}".format(source))
    
    assert sum([1 for symbol in source if symbol.startswith(Symbols.DIS.value)]) == sum([1 for symbol in target if symbol.startswith(Symbols.DIS.value)]), "Number of dispersion events must be equal between source and target sequence"


def remove_file(file):
    if os.path.exists(file):
        os.remove(file)


def reset_file(filename):
    with open(filename, "w") as f_reset:
        f_reset.truncate()


def generate_seq(length):
    base_map = {1: "A", 2: "T", 3: "G", 4: "C"}
    return ''.join([base_map[random.randint(1, 4)] for _ in range(length)])


def percent_N(seq):
    return 0 if len(seq) == 0 else seq.count('N') / len(seq)


def complement(seq):
    output = ""
    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for base in seq.upper():
        if base in base_complements:
            output += base_complements[base]
        else:
            output += base

    return output


def reformat_seq(transformation):
    # transformation: tuple, user inputted source and target
    # if dispersion events exist in transformation, tag on unique ids to make them distinct as they all are "_"
    unique_transform = []
    unique_id = 1
    for component in transformation:
        if component != Symbols.DIS.value and component != Symbols.DUP.value and component != Symbols.DIV.value:
            unique_transform.append(component)
        elif component == Symbols.DUP.value:  # duplication event case, need to group together symbol and duplication marking
            unique_transform[-1] += Symbols.DUP.value
        elif component == Symbols.DIV.value:  # divergence event case, want to keep track of interval needing modification
            unique_transform[-1] += Symbols.DIV.value
        else:  # dispersion event case, component = dispersion
            unique_transform.append(component + str(unique_id))
            unique_id += 1
    return tuple(unique_transform)


def divergence(seq, divergence_prob=None):
    # apply random base flips to input sequence
    # probability of changing the base, chosen from U(0.5,1.0) or given from user
    p = random.uniform(0.5, 1.0) if divergence_prob is None else float(divergence_prob)
    return ''.join([b if random.random() > p else random.choice(list({"A", "C", "T", "G"} - {b})) for b in seq.upper()])


def get_sv_config_identifier(sv_config):
    return sv_config['variant_set_id']


class RegionsOfInterest:
    """
    container class for the optionally-provided genome context elements that will be used in SV placement
    """
    def __init__(self):
        self.roi_dict = defaultdict(list)

    def parse_bed_file(self, bed_fname, allow_chroms=None):
        """
        extracts bed file intervals to be considered for randomized event overlap and adds
        (chr, start, end) tuples representing the intervals to a dict keyed on elt type
        - allow_chroms: optional list of allowed chromosomes (filtering out all entries with chrom not in list)
        """
        with open(bed_fname, 'r') as bed_file:
            for line in bed_file:
                if line.startswith('#') or line.isspace():
                    continue
                fields = line.strip().split()
                assert len(fields) >= 3, "Unexpected number of fields in BED: %s" % line
                chr_name, start, end, elt_type = fields[:4]
                if allow_chroms and chr_name not in allow_chroms:
                    continue
                self.roi_dict[elt_type].append((chr_name, int(start), int(end)))
        # shuffle the lists of each element type for random selection
        for k in self.roi_dict.keys():
            random.shuffle(self.roi_dict[k])

    @staticmethod
    def midpoint(start, end):
        return (start + end) // 2

    @staticmethod
    def get_intrvl_len(chr, st, end):
        return int(end) - int(st)


class BlacklistRegions(RegionsOfInterest):
    def __init__(self, config):
        super(BlacklistRegions, self).__init__()

        if isinstance(config['blacklist_regions'], list):
            for path in config['blacklist_regions']:
                self.parse_file(path)
        else:
            self.parse_file(config['blacklist_regions'])

    def parse_file(self, path):
        if path[-4:] == '.vcf':
            self.parse_vcf_file(path)
        elif path[-4:] == '.bed':
            self.parse_bed_file(path)
        else:
            raise Exception('Only .bed and .vcf files may be given in blacklist_regions')

    def parse_vcf_file(self, vcf_fname):
        for rec in VariantFile(vcf_fname):
            self.roi_dict[rec.id].append((rec.chrom, rec.start, rec.stop))

    def get_roi(self, sv_config):
        """
        For a given sv_config entry, return all applicable
        blacklist_regions that will have to be avoided at SV placement
        """
        if sv_config['blacklist_region_type'].lower() == 'all':
            return list(itertools.chain(*self.roi_dict.values()))
        else:
            # if type given, get regions whose types are prefixed by specified string
            return list(itertools.chain(*[v for k, v in self.roi_dict.items() if
                        k.startswith(sv_config['blacklist_region_type'])]))


class OverlapRegions(RegionsOfInterest):
    def __init__(self, config, allow_chroms=None):
        # elements will be stored in a dict keyed on element name
        super().__init__()
        self.allow_chroms = allow_chroms
        self.svtype_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_partial_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_alu_mediated_counts = defaultdict(int)
        self.get_num_overlap_counts(config)

        if isinstance(config['overlap_regions'], list):
            for bed_path in config['overlap_regions']:
                self.parse_bed_file(bed_path, allow_chroms=self.allow_chroms)
        else:
            self.parse_bed_file(config['overlap_regions'], allow_chroms=self.allow_chroms)

        # optional list of Alu pairs to be used as flanking Alus for DELs/DUPs (to simulate Alu-mediated CNVs)
        # --> dict of form {SV_identifier: [(chrom_1, a_1, b_1), ..., (chrom_n, a_n, b_n)]}
        self.alu_pairs = self.populate_alu_pairs(config['variant_sets'])

    def get_num_overlap_counts(self, config):
        # populate nested dict of the form {sv type: {element type: num_overlap}}
        for sv in config['variant_sets']:
            if 'overlap_type' not in sv:
                if ('overlap_component' in sv and sv['overlap_component'].lower() == 'target') or \
                        ('overlap_region_type' in sv and isinstance(sv['overlap_region_type'], list)):
                    # if placing a target locus in a ROI, overlap type can be omitted and will be treated
                    # as 'exact' for the purpose of interval selection (since specific target locus is
                    # chosen from the selected ROI at SV initialization)
                    sv['overlap_type'] = 'exact'
                else:
                    continue
            elif 'overlap_component' in sv and sv['overlap_component'].lower() == 'target':
                raise Exception('CONFIG ERROR: if \'overlap_component\': \'target\', must omit \'overlap_type\'')
            sv_config_key = get_sv_config_identifier(sv)
            # building separate counts dictionaries for complete overlaps, partial overlaps, and alu-mediated intervals
            for overlap_type, count_dict in [('exact', self.svtype_overlap_counts), ('partial', self.svtype_partial_overlap_counts)]:
                if sv['overlap_type'] == overlap_type:
                    # set overlap type label based on specified overlap_region_type
                    if 'overlap_region_type' in sv:
                        if isinstance(sv['overlap_region_type'], str):
                            count_dict[sv_config_key][sv['overlap_region_type']] += sv['number']
                        elif isinstance(sv['overlap_region_type'], list):
                            for roi_type in sv['overlap_region_type']:
                                if roi_type is not None:
                                    count_dict[sv_config_key][roi_type] += sv['number']
                        else:
                            raise Exception('\'overlap_region_type\' must be of type string, or a list of strings')
                    else:
                        count_dict[sv_config_key]['ALL'] += sv['number']
            if sv['overlap_type'] == 'flanked':
                assert 'overlap_region_type' in sv and sv['overlap_region_type'] == 'Alu', \
                    'For \'overlap_type\': \'flanked\', must specify \'overlap_region_type\': \'Alu\''
                self.svtype_alu_mediated_counts[sv_config_key] += sv['number']

    def get_overlap_intervals(self, sv_config_id, sv_config, partial_overlap):
        if len(self.roi_dict) == 0:
            return None
        frags = FormatterIO.get_grammar(sv_config)
        counts_dict = self.svtype_partial_overlap_counts if partial_overlap else self.svtype_overlap_counts
        elt_types = list(counts_dict[sv_config_id].keys())
        for elt_type in elt_types:
            # if this elt_type doesn't have any matching elts in the element dict, remove it and move on to the next elt
            while elt_type != 'ALL' and not any([elt_type in ovlp_type for ovlp_type in self.roi_dict.keys()]):
                del counts_dict[sv_config_id][elt_type]
                if len(counts_dict[sv_config_id]) == 0:
                    del counts_dict[sv_config_id]
                    break
        if 'overlap_component' not in sv_config or sv_config['overlap_component'] == 'source':
            # 'overlap_component' will be left out if either SV is simple or a list is given for 'overlap_region_type'
            if 'overlap_region_type' in sv_config and isinstance(sv_config['overlap_region_type'], list):
                # if overlap_region_type given as list of elt_types in positions corresponding to their matching fragments
                # pull those lists apart so they can be traversed together at interval retrieval time
                frag_idx, elt_types = list(zip(*[(i, v) for i, v in enumerate(sv_config['overlap_region_type']) if v is not None]))
            else:
                frag_idx = [0]
        elif sv_config['overlap_component'] == 'target':
            frag_idx = [1]
        elif sv_config['overlap_component'] == 'rand':
            frag_idx = [random.randint(0, len(sv_config['length_ranges']) - 1)]
        elif sv_config['overlap_component'] == 'full_sv':
            frag_idx = list(range(len(sv_config['length_ranges'])))
        else:
            raise Exception('Invalid \'overlap_component\' argument: %s' % sv_config['overlap_component'])
        repeat_elts = {}
        if 'overlap_component' in sv_config and sv_config['overlap_component'] == 'full_sv':
            elt_type = elt_types[0]
            # sum over all involved frags in case overlap_component == 'full_sv'
            minsize = sum(sv_config['length_ranges'][i][0] for i in frag_idx)
            maxsize = sum(sv_config['length_ranges'][i][1] for i in frag_idx)
            repeat_elt = self.get_roi(sv_config_id=sv_config_id,
                                      minsize=minsize, maxsize=maxsize,
                                      elt_type=(None if elt_type == 'ALL' else elt_type),
                                      partial_overlap=partial_overlap)
            repeat_elts['full_sv'] = repeat_elt
        else:
            for idx, f_id in enumerate(frag_idx):
                if frags[f_id].startswith(Symbols.DIS.value):
                    minsize, maxsize = None, None
                else:
                    minsize = sv_config['length_ranges'][f_id][0]
                    maxsize = sv_config['length_ranges'][f_id][1]
                elt_type = elt_types[0] if len(elt_types) == 1 else elt_types[idx]
                repeat_elt = self.get_roi(sv_config_id=sv_config_id,
                                          minsize=minsize, maxsize=maxsize,
                                          elt_type=(None if elt_type == 'ALL' else elt_type),
                                          partial_overlap=partial_overlap)
                repeat_elts[frags[f_id]] = repeat_elt
        return repeat_elts

    def populate_alu_pairs(self, svs_config):
        if not any(('overlap_type', 'flanked') in d.items() and ('overlap_region_type', 'Alu') in d.items() for d in svs_config):
            return None
        # construct dict of form {chrom: [(start, end)]} Alu interval tuples
        alu_intervals = defaultdict(list)
        for elt_type, elt_list in self.roi_dict.items():
            if 'Alu' in elt_type:
                for chrom, start, end in elt_list:
                    alu_intervals[chrom].append((int(start), int(end)))
        for chrom in alu_intervals:
            random.shuffle(alu_intervals[chrom])
        alu_pairs_dict = {}
        # construct the lists of viable flanking pairs for each SV config with a nonzero num_alu_mediated specified
        for sv_config in svs_config:
            sv_config_id = get_sv_config_identifier(sv_config)
            if ('overlap_type', 'flanked') in sv_config.items() and ('overlap_region_type', 'Alu') in sv_config.items():
                alu_intervals_copy = copy.deepcopy(alu_intervals)
                alu_pairs = []
                sv_min, sv_max = sv_config['length_ranges'][0]  # sv_config['min_length'], sv_config['max_length']
                for _ in range(self.svtype_alu_mediated_counts[sv_config_id]):
                    viable_matches = []
                    while len(viable_matches) == 0:
                        if len(alu_intervals_copy) == 0 or not any(len(v) > 1 for v in alu_intervals_copy.values()):
                            break  # <-- no Alus remaining to choose from
                        # choose a random starting (left) Alu
                        rand_chrom = random.choice([chrom for chrom, alus in alu_intervals_copy.items() if len(alus) > 1])
                        alu1 = alu_intervals_copy[rand_chrom].pop()
                        # collect all the alus on the same chromosome that are within an appropriate dist
                        viable_matches = [ivl for ivl in alu_intervals_copy[rand_chrom] if sv_min <= np.abs(self.midpoint(*ivl) - self.midpoint(*alu1)) <= sv_max]
                    if len(viable_matches) > 0:
                        alu2 = random.choice(viable_matches)
                        alu_intervals_copy[rand_chrom].remove(alu2)
                        if len(alu_intervals_copy[rand_chrom]) == 0:
                            del alu_intervals_copy[rand_chrom]
                        alu_pairs.append((rand_chrom, alu1, alu2) if alu1[0] < alu2[0] else (rand_chrom, alu2, alu1))
                        # delete the used Alus from the alu_intervals dict and from the general roi_dict
                        for alu in [alu1, alu2]:
                            alu_intervals[rand_chrom].remove(alu)
                            self.remove_alu_from_overlap_dict(rand_chrom, *alu)
                if len(alu_pairs) > 0:
                    random.shuffle(alu_pairs)
                    alu_pairs_dict[sv_config_id] = alu_pairs
        return alu_pairs_dict

    def get_alu_mediated_interval(self, sv_config_id):
        if len(self.alu_pairs[sv_config_id]) == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
            return None, None, None
        ivl_chrom, left_alu, right_alu = self.alu_pairs[sv_config_id].pop()
        ivl_start, ivl_end = (left_alu[1] + left_alu[0]) // 2, (right_alu[1] + right_alu[0]) // 2
        self.svtype_alu_mediated_counts[sv_config_id] -= 1
        if self.svtype_alu_mediated_counts[sv_config_id] == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
        # returning overlap element and type in dictionaries to match output format of get_overlap_intervals()
        return {'A': (ivl_chrom, ivl_start, ivl_end, 'ALU_MEDIATED')}

    def remove_alu_from_overlap_dict(self, chrom, start, end):
        # helper method to remove an Alu element from the roi_dict if it's been chosen for an alu-mediated interval
        for k, v in self.roi_dict.items():
            if k.startswith('Alu') and (chrom, start, end) in v:
                self.roi_dict[k].remove((chrom, start, end))
        self.roi_dict = {k: v for k, v in self.roi_dict.items() if len(v) > 0}

    @staticmethod
    def get_partially_overlapping_interval(elt_chrom, elt_start, elt_stop, sv_min, sv_max):
        # choose sv size, draw a random position in element interval, randomly set that to the start or end of SV
        def draw_from_unif(a, b):
            return a if a == b else np.random.randint(a, b)
        sv_size = draw_from_unif(sv_min, sv_max)
        interval_left = draw_from_unif(max(elt_start - sv_size + 1, 0), max(elt_stop - sv_size, 0)) if random.randint(0, 1)\
            else draw_from_unif(elt_start, max(elt_stop - sv_size - 1, elt_start))  # <- two possible distributions for the left bound of the SV
        return elt_chrom, interval_left, interval_left + sv_size

    def decrement_counts(self, sv_config_id, input_elt_type, partial_overlap):
        # decrement the relevant counts dict based on element extracted in get_repeat_elt()
        counts_dict = self.svtype_partial_overlap_counts if partial_overlap else self.svtype_overlap_counts
        counts_dict[sv_config_id][input_elt_type] -= 1
        if counts_dict[sv_config_id][input_elt_type] == 0:
            del counts_dict[sv_config_id][input_elt_type]
            if len(counts_dict[sv_config_id]) == 0:
                del counts_dict[sv_config_id]

    def get_roi(self, sv_config_id, minsize, maxsize, elt_type=None, partial_overlap=False):
        input_elt_type = 'ALL' if elt_type is None else elt_type
        if elt_type in self.roi_dict.keys():  # <- elt_type given, and elements of matching type are in the dict
            rand_elt = next((elt for elt in self.roi_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None) \
                        if not partial_overlap and minsize is not None and maxsize is not None \
                        else next(elt for elt in self.roi_dict[elt_type])
        elif len(self.roi_dict.keys()) > 0:  # <- element dict is non-empty; elt_type not given, or is given as a prefix
            if elt_type is None:
                elt_type_mapping = {elt: elt_type for (elt_type, elt_list) in self.roi_dict.items() for elt in elt_list}
                # elt_type-agnostic branch: draw random element from combined list of elements across all types
                rand_elt = next((elt for elt in itertools.chain.from_iterable(self.roi_dict.values()) if
                                 minsize <= self.get_intrvl_len(*elt) <= maxsize), None) \
                    if not partial_overlap and minsize is not None and maxsize is not None \
                    else next((elt for elt in itertools.chain.from_iterable(self.roi_dict.values())), None)
            else:
                # elt_type prefix branch: choose an element whose name begins with elt_type
                prefix_matches = [element_type for element_type in self.roi_dict.keys() if element_type.startswith(elt_type)]
                if len(prefix_matches) > 0:
                    random.shuffle(prefix_matches)
                    elt_type = prefix_matches.pop()
                    while len(prefix_matches) > 0 and len([elt for elt in self.roi_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize]) == 0:
                        elt_type = prefix_matches.pop()
                    rand_elt = next((elt for elt in self.roi_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None) \
                        if not partial_overlap and minsize is not None and maxsize is not None \
                        else next(elt for elt in self.roi_dict[elt_type])
                else:
                    rand_elt = None
        else:  # <- elt_type is NOT given and there AREN'T elements in the dict
            rand_elt = None
        if rand_elt is not None:
            if elt_type is None:
                elt_type = elt_type_mapping[rand_elt]
            self.roi_dict[elt_type].remove(rand_elt)
            if len(self.roi_dict[elt_type]) == 0:
                del self.roi_dict[elt_type]
            self.decrement_counts(sv_config_id, input_elt_type, partial_overlap)
        elif elt_type in self.roi_dict.keys():  # <- case in which elements of valid type but not size were available
            self.decrement_counts(sv_config_id, input_elt_type, partial_overlap)
        if partial_overlap and rand_elt is not None:
            rand_elt = self.get_partially_overlapping_interval(*rand_elt, minsize, maxsize)
        return (*rand_elt, elt_type) if rand_elt is not None else None
