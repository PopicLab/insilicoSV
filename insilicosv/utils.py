import os
import random
import itertools
import copy  # <- utils for making copies of python objects
import argparse
import numpy as np
from insilicosv.constants import *
from collections import defaultdict
from pysam import VariantFile


class NestedDict(defaultdict):
    def __call__(self):
        return NestedDict(self.default_factory)


def get_original_base_pos(genome_vcf, query_loc, query_chrom):
    # calculates the input reference position of the base that corresponds
    # to the one at position query_loc in the output reference (i.e., adjusts the
    # input ref position by the cumulative change in ref length caused by
    # the variants from [0, query_loc)
    # **for now: just admits DEL, DUP, INV
    vcf = VariantFile(genome_vcf)
    # need to make sure vcf recs sorted by position and separated by chrom
    vcf_recs = defaultdict(list)
    for rec in vcf.fetch():
        vcf_recs[rec.chrom].append((rec.id, rec.start, rec.stop))
    for chrom in vcf_recs.keys():
        vcf_recs[chrom].sort(key=lambda x: x[1])
    # total the cumulative change in ref length from all the SVs *ending* before p
    p = query_loc
    i = 0
    while i < len(vcf_recs[query_chrom]) and vcf_recs[query_chrom][i][1] <= p:
        rec_id, rec_start, rec_stop = vcf_recs[query_chrom][i]
        rec_len = rec_stop - rec_start
        if rec_id == 'DEL' and rec_start <= p:
            p += rec_len
        if rec_id == 'DUP' and rec_stop <= p:
            p -= rec_len
        i += 1
    # check if p is internal to an INV; edge case not covered by the above loop
    surrounding_inv = next((rec for rec in vcf_recs[query_chrom] if rec[1] <= p < rec[2] and rec[0] == 'INV'), None)
    if surrounding_inv is not None:
        sv_id, sv_start, sv_end = surrounding_inv
        p += (sv_end - sv_start) - ((p - sv_start) * 2 + 1)
    return p


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


def divergence(seq, divergence_prob=None):
    # apply random base flips to input sequence
    # probability of changing the base, chosen from U(0.5,1.0) or given from user
    p = random.uniform(0.5, 1.0) if divergence_prob is None else float(divergence_prob)
    return ''.join([b if random.random() > p else random.choice(list({"A", "C", "T", "G"} - {b})) for b in seq.upper()])


def get_sv_config_identifier(sv_config):
    # helper function to yield a uniquely identifying string for a given SV config entry
    # --> usage: needing to index into overlap events dicts to access info specific to a set of SVs simulated with a given set of params
    min_str = str(sv_config['min_length']) if isinstance(sv_config['min_length'], int) else '_'.join(map(str, sv_config['min_length']))
    max_str = str(sv_config['max_length']) if isinstance(sv_config['max_length'], int) else '_'.join(map(str, sv_config['max_length']))
    return sv_config['type'].value + '_' + min_str + '_' + max_str


class OverlapEvents:
    # container object for the optionally-provided genome context elements that will be used in SV placement
    def __init__(self, config, allow_chroms=None):
        # elements will be stored in a dict keyed on element name
        self.overlap_events_dict = defaultdict(list)
        self.allow_chroms = allow_chroms
        self.allow_types = None if 'allow_types' not in config['overlap_events'].keys() else config['overlap_events']['allow_types']
        if isinstance(self.allow_types, str):
            self.allow_types = [self.allow_types]
        self.svtype_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_partial_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_alu_mediated_counts = defaultdict(NestedDict(int))
        self.get_num_overlap_counts(config)

        if type(config['overlap_events']['bed']) is list:
            for bed_path in config['overlap_events']['bed']:
                self.parse_bed_file(bed_path, allow_chroms=self.allow_chroms, allow_types=self.allow_types)
        else:
            self.parse_bed_file(config['overlap_events']['bed'], allow_chroms=self.allow_chroms, allow_types=self.allow_types)

        # optional list of Alu pairs to be used as flanking Alus for DELs/DUPs (to simulate Alu-mediated CNVs)
        # --> dict of form {SV_identifier: [(chrom_1, a_1, b_1), ..., (chrom_n, a_n, b_n)]}
        self.alu_pairs = self.populate_alu_pairs(config['SVs'])

    def get_num_overlap_counts(self, config):
        # populate nested dict of the form {sv type: {element type: num_overlap}}
        for sv in config['SVs']:
            sv_config_key = get_sv_config_identifier(sv)
            # building separate counts dictionaries for complete overlaps, partial overlaps, and alu-mediated intervals
            for overlap_count, count_dict in [('num_overlap', self.svtype_overlap_counts), ('num_partial_overlap', self.svtype_partial_overlap_counts)]:
                if overlap_count in sv.keys():
                    if isinstance(sv[overlap_count], int) and sv[overlap_count] > 0:
                        # set overlap type label based on specified allow_types
                        if self.allow_types is not None and len(self.allow_types) == 1:
                            count_dict[sv_config_key][self.allow_types[0]] = sv[overlap_count]
                        else:
                            count_dict[sv_config_key]['ALL'] = sv[overlap_count]
                    else:
                        # define correspondence with the element types as they were given in the config file
                        for i in range(len(sv[overlap_count])):
                            if sv[overlap_count][i] > 0:  # <- only counting entries in list > 0
                                count_dict[sv_config_key][self.allow_types[i]] = sv[overlap_count][i]
            if 'num_alu_mediated' in sv.keys():
                self.svtype_alu_mediated_counts[sv_config_key] = sv['num_alu_mediated']

    def parse_bed_file(self, bed_fname, allow_chroms=None, allow_types=None):
        """
        reads bed file (intended use: processing repeatmasker elements to be considered for randomized
        event overlap) and adds (chr, start, end) tuples representing the intervals to a dict keyed on elt type
        - allow_chroms: optional list of allowed chromosomes (filtering out all entries with chrom not in list)
        - allow_types: optional list of allowed types to be included in the simulation
        """
        with open(bed_fname, 'r') as bed_file:
            for line in bed_file:
                if line.startswith('#') or line.isspace():
                    continue
                fields = line.strip().split()
                assert len(fields) >= 3, "Unexpected number of fields in BED: %s" % line
                chr_name, start, end, elt_type = fields[:4]
                # need to allow for allow_types to include prefixes of element names (always including Alu as an allowed
                # prefix in case any of the SVs will need to be placed at Alu-mediated intervals)
                if (allow_chroms and chr_name not in allow_chroms) or (allow_types and not any([elt_name in elt_type for elt_name in allow_types + ['Alu']])):
                    continue
                self.overlap_events_dict[elt_type].append((chr_name, int(start), int(end)))
        # shuffle the lists of each element type for random selection
        for k in self.overlap_events_dict.keys():
            random.shuffle(self.overlap_events_dict[k])

    def get_single_element_interval(self, sv_config_id, sv_config, partial_overlap):
        if len(self.overlap_events_dict) == 0:
            return None, None, None
        counts_dict = self.svtype_partial_overlap_counts if partial_overlap else self.svtype_overlap_counts
        elt_type = list(counts_dict[sv_config_id].keys())[0]
        # if this elt_type doesn't have any matching elts in the element dict, remove it and move on to the next elt
        while elt_type != 'ALL' and not any([elt_type in ovlp_type for ovlp_type in self.overlap_events_dict.keys()]):
            del counts_dict[sv_config_id][elt_type]
            if len(counts_dict[sv_config_id]) == 0:
                del counts_dict[sv_config_id]
                break
            elt_type = list(counts_dict[sv_config_id].keys())[0]
        # if num_overlaps was given as a single number rather than an element-specific list, draw a random element for overlap
        repeat_elt, retrieved_type = self.__getitem__(sv_config_id=sv_config_id,
                                                      minsize=sv_config['length_ranges'][0][0],
                                                      maxsize=sv_config['length_ranges'][0][1],
                                                      elt_type=(None if elt_type == 'ALL' else elt_type),
                                                      partial_overlap=partial_overlap)
        return repeat_elt, retrieved_type, elt_type

    def populate_alu_pairs(self, svs_config):
        if not any('num_alu_mediated' in d.keys() for d in svs_config):
            return None
        # construct dict of form {chrom: [(start, end)]} Alu interval tuples
        alu_intervals = defaultdict(list)
        for elt_type, elt_list in self.overlap_events_dict.items():
            if 'Alu' in elt_type:
                for chrom, start, end in elt_list:
                    alu_intervals[chrom].append((int(start), int(end)))
        for chrom in alu_intervals.keys():
            random.shuffle(alu_intervals[chrom])
        alu_pairs_dict = {}
        # construct the lists of viable flanking pairs for each SV config with a nonzero num_alu_mediated specified
        for sv_config in svs_config:
            sv_config_id = get_sv_config_identifier(sv_config)
            if 'num_alu_mediated' in sv_config.keys():
                alu_intervals_copy = copy.deepcopy(alu_intervals)
                alu_pairs = []
                sv_min, sv_max = sv_config['min_length'], sv_config['max_length']
                for _ in range(self.svtype_alu_mediated_counts[sv_config_id]):
                    viable_matches = []
                    while len(viable_matches) == 0:
                        if len(alu_intervals_copy.keys()) == 0 or not any(len(v) > 1 for v in alu_intervals_copy.values()):
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
                        # delete the used Alus from the alu_intervals dict and from the general overlap_events_dict
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
            return None, None
        ivl_chrom, left_alu, right_alu = self.alu_pairs[sv_config_id].pop()
        ivl_start, ivl_end = (left_alu[1] + left_alu[0]) // 2, (right_alu[1] + right_alu[0]) // 2
        self.svtype_alu_mediated_counts[sv_config_id] -= 1
        if self.svtype_alu_mediated_counts[sv_config_id] == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
        return (ivl_chrom, ivl_start, ivl_end), 'ALU_MEDIATED'  # <- ALU_MEDIATED to be used as retrieved_type

    def remove_alu_from_overlap_dict(self, chrom, start, end):
        # helper method to remove an Alu element from the overlap_events_dict if it's been chosen for an alu-mediated interval
        for k, v in self.overlap_events_dict.items():
            if k.startswith('Alu') and (chrom, start, end) in v:
                self.overlap_events_dict[k].remove((chrom, start, end))
        self.overlap_events_dict = {k: v for k, v in self.overlap_events_dict.items() if len(v) > 0}

    @staticmethod
    def midpoint(start, end):
        return (start + end) // 2

    @staticmethod
    def get_intrvl_len(chr, st, end):
        return int(end) - int(st)

    def elt_type_is_allowed(self, elt_type):
        # helper to check if an element type is allowed by the allow_types list
        return True if self.allow_types is None else any(allow_type.startswith(elt_type) for allow_type in self.allow_types)

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
        # decrement the relevant counts dict based on element extracted in __getitem__()
        counts_dict = self.svtype_partial_overlap_counts if partial_overlap else self.svtype_overlap_counts
        counts_dict[sv_config_id][input_elt_type] -= 1
        if counts_dict[sv_config_id][input_elt_type] == 0:
            del counts_dict[sv_config_id][input_elt_type]
            if len(counts_dict[sv_config_id].keys()) == 0:
                del counts_dict[sv_config_id]

    def __getitem__(self, sv_config_id, minsize, maxsize, elt_type=None, partial_overlap=False):
        input_elt_type = 'ALL' if elt_type is None else elt_type
        if elt_type in self.overlap_events_dict.keys():  # <- elt_type given, and elements of matching type are in the dict
            rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None) \
                        if not partial_overlap else next(elt for elt in self.overlap_events_dict[elt_type])
        elif len(self.overlap_events_dict.keys()) > 0:  # <- element dict is non-empty; elt_type not given, or is given as a prefix
            if elt_type is None:
                elt_type_mapping = {elt: elt_type for (elt_type, elt_list) in self.overlap_events_dict.items() for elt in elt_list}
                # elt_type-agnostic branch: draw random element from combined list of elements across all types
                # and if allow_types is specified, only grab agreeing elements
                rand_elt = next((elt for elt in itertools.chain.from_iterable(self.overlap_events_dict.values()) if
                                 (minsize <= self.get_intrvl_len(*elt) <= maxsize or partial_overlap) and
                                 (self.elt_type_is_allowed(elt_type_mapping[elt]))), None)
            else:
                # elt_type prefix branch: choose an element whose name begins with elt_type
                prefix_matches = [element_type for element_type in self.overlap_events_dict.keys() if element_type.startswith(elt_type)]
                if len(prefix_matches) > 0:
                    random.shuffle(prefix_matches)
                    elt_type = prefix_matches.pop()
                    while len(prefix_matches) > 0 and len([elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize]) == 0:
                        elt_type = prefix_matches.pop()
                    rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None) \
                        if not partial_overlap else next(elt for elt in self.overlap_events_dict[elt_type])
                else:
                    rand_elt = None
        else:  # <- elt_type is NOT given and there AREN'T elements in the dict
            rand_elt = None
        if rand_elt is not None:
            if elt_type is None:
                elt_type = elt_type_mapping[rand_elt]
            self.overlap_events_dict[elt_type].remove(rand_elt)
            if len(self.overlap_events_dict[elt_type]) == 0:
                del self.overlap_events_dict[elt_type]
            self.decrement_counts(sv_config_id, input_elt_type, partial_overlap)
        elif elt_type in self.overlap_events_dict.keys():  # <- case in which elements of valid type but not size were available
            self.decrement_counts(sv_config_id, input_elt_type, partial_overlap)
        if partial_overlap and rand_elt is not None:
            rand_elt = self.get_partially_overlapping_interval(*rand_elt, minsize, maxsize)
        return rand_elt, elt_type


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Special use case for utils.py: Commandline access to get_original_base_pos()')
    parser.add_argument('--genome_vcf', type=str, help='Path to vcf describing simulated SVs populating the synthetic genome')
    parser.add_argument('--query_loc', type=int, help='Query genome locus (given with respect to synthetic reference)')
    parser.add_argument('--query_chrom', type=int, help='Corresponding chromosome of the query locus')
    args = parser.parse_args()

    print('Input reference position corresponding to mutated reference position (%s) %d: %d' %
          (args.query_chrom, args.query_loc, get_original_base_pos(args.genome_vcf, args.query_loc, args.query_chrom)))

