import os
import random
import itertools
import copy  # <- utils for making copies of python objects
import argparse
import numpy as np
from constants import *
from collections import defaultdict
from pysam import VariantFile


# NestedDict helper class
class NestedDict(defaultdict):
    """Nested defaultdict"""

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
    # print(f'\ninput query locus: p = {p}')
    while i < len(vcf_recs[query_chrom]) and vcf_recs[query_chrom][i][1] <= p:
        rec_id, rec_start, rec_stop = vcf_recs[query_chrom][i]
        rec_len = rec_stop - rec_start
        if rec_id == 'DEL' and rec_start <= p:
            p += rec_len
            # print(f'DEL found, adding {rec_len}; p = {p}')
        if rec_id == 'DUP' and rec_stop <= p:
            p -= rec_len
            # print(f'DUP found, subtracting {rec_len}; p = {p}')
        i += 1
    # print(f'accounted for prior SVs:\tp = {p}')
    # ... then check if p is internal to an INV; edge case not covered by the above loop
    surrounding_inv = next((rec for rec in vcf_recs[query_chrom] if rec[1] <= p < rec[2] and rec[0] == 'INV'), None)
    if surrounding_inv is not None:
        sv_id, sv_start, sv_end = surrounding_inv
        p += (sv_end - sv_start) - ((p - sv_start) * 2 + 1)
    # print(f'input: {query_loc}\toutput: {p}')
    return p


def is_overlapping(event_ranges, addition, called_from_helper=False, strictly_partial=False):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored
    # --> 'strictly_partial' a toggle to return True if the interval overlaps a stored event, without being equal to it
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
    # removes file if it exists
    if os.path.exists(file):
        os.remove(file)

def reset_file(filename):
    #print("Overwritting File {}...".format(filename))
    with open(filename, "w") as f_reset:
        f_reset.truncate()

def generate_seq(length):
    # helper function for insertions
    # generates random sequence of bases of given length
    base_map = {1: "A", 2: "T", 3: "G", 4: "C"}
    return ''.join([base_map[random.randint(1, 4)] for _ in range(length)])

def percent_N(seq):
    return 0 if len(seq) == 0 else seq.count('N') / len(seq)

def complement(seq):
    # seq: str
    output = ""
    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for base in seq.upper():
        if base in base_complements:
            output += base_complements[base]
        else:
            output += base
            # raise ValueError("Unknown base \'{}\' detected in reference, complement of base not taken".format(base))

    return output

def divergence(seq):
    # function to create slightly mutated version of an input sequence
    # --> p given as the probability of changing the base, chosen from U(0.5,1.0)
    p = random.uniform(0.5,1.0)
    return ''.join([b if random.random() > p else random.choice(list({"A", "C", "T", "G"} - {b})) for b in seq.upper()])


def get_sv_config_identifier(sv_config):
    # helper function to yield a uniquely identifying string
    # for a given SV config entry
    # --> usage: needing to index into overlap events dicts to access info specific to a set of SVs simulated with a given set of params
    min_str = str(sv_config['min_length']) if isinstance(sv_config['min_length'], int) else '_'.join(map(str, sv_config['min_length']))
    max_str = str(sv_config['max_length']) if isinstance(sv_config['max_length'], int) else '_'.join(map(str, sv_config['max_length']))
    return sv_config['type'].value + '_' + min_str + '_' + max_str


# container object for the optionally-provided genome context elements that will be used in SV placement
class OverlapEvents:
    def __init__(self, config, allow_chroms=None):
        # elements will be stored in a dict keyed on element name
        self.overlap_events_dict = defaultdict(list)
        self.allow_chroms = allow_chroms
        # --> allow_types will optionally be given in the config
        self.allow_types = None if 'allow_types' not in config['overlap_events'].keys() else config['overlap_events']['allow_types']
        # if a single allow_type is given, wrap it in a list for agreement with downstream logic
        if isinstance(self.allow_types, str):
            self.allow_types = [self.allow_types]
        # dicts containing the counts of overlapping SVs and alu_mediated SVs for each SV type and known element type
        self.svtype_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_partial_overlap_counts = defaultdict(NestedDict(int))
        self.svtype_alu_mediated_counts = defaultdict(NestedDict(int))
        self.get_num_overlap_counts(config)
        # # ---- debug ----
        # print(f'overlap_counts: {self.svtype_overlap_counts}')
        # print(f'partial overlap_counts: {self.svtype_partial_overlap_counts}')
        # print(f'alu_mediated_counts: {self.svtype_alu_mediated_counts}')
        # # ---------------

        # overlap_events can either be given as a single .bed file or multiple
        if type(config['overlap_events']['bed']) is list:
            for bed_path in config['overlap_events']['bed']:
                self.parse_bed_file(bed_path, allow_chroms=self.allow_chroms, allow_types=self.allow_types)
        else:
            self.parse_bed_file(config['overlap_events']['bed'], allow_chroms=self.allow_chroms, allow_types=self.allow_types)
        # # ---- debug ----
        # print(f'overlap_events: {self.overlap_events_dict}')
        # # ---------------

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
                        # if num_overlap given as a singleton integer, label type as 'ALL' if zero or multiple allow_types
                        # are given, otherwise set type to the single allow_type provided
                        if self.allow_types is not None and len(self.allow_types) == 1:
                            count_dict[sv_config_key][self.allow_types[0]] = sv[overlap_count]
                        else:
                            count_dict[sv_config_key]['ALL'] = sv[overlap_count]
                    else:
                        # set up correspondence with the element types as they were given in the config file
                        for i in range(len(sv[overlap_count])):
                            if sv[overlap_count][i] > 0:  # <- only counting entries in list > 0
                                count_dict[sv_config_key][self.allow_types[i]] = sv[overlap_count][i]
            if 'num_alu_mediated' in sv.keys():
                self.svtype_alu_mediated_counts[sv_config_key] = sv['num_alu_mediated']

    def parse_bed_file(self, bed_fname, allow_chroms=None, allow_types=None):
        """
        reads bed file (intended use: processing repeatmasker elements to be considered for randomized
        event overlap) and adds (chr, start, end) tuples representing the intervals to a dict keyed on elt type
        --> logic taken from bed_iter() and parse_bed_line() in cue (seq/io.py)
        - allow_chroms: optional list of allowed chromosomes (filtering out all entries with chrom not in list)
        - allow_types: optional list of allowed types to be included in the simulation
        """
        with open(bed_fname, 'r') as bed_file:
            for line in bed_file:
                if line.startswith('#') or line.isspace():
                    continue
                # parse line
                fields = line.strip().split()
                assert len(fields) >= 3, "Unexpected number of fields in BED: %s" % line
                chr_name, start, end, elt_type = fields[:4]
                # need to allow for allow_types to include prefixes of element names
                # --> also need to always include Alu as an allowed prefix in case any of the SVs will need to be placed at Alu-mediated intervals
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
        # helper method to find a whole or partial overlap interval
        elt_type = list(counts_dict[sv_config_id].keys())[0]
        # # ---- debug ----
        # print(f'\nget_single_interval() called with:\n--> sv_config_id: {sv_config_id}\n--> counts_dict: {counts_dict}\n--> elt_type: {elt_type}'
        #       f'\n--> candidate_elts = {self.overlap_events_dict}')
        # # ---------------
        # if this elt_type doesn't have any matching elements in the element dict (for matching determined
        # either with elt_type giving the full element name or just a prefix – e.g., 'L1'), remove it
        # and move on to the next elt type
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
        # alu_pairs only to be populated if any SV config entries specify a value for 'num_alu_mediated'
        if not any('num_alu_mediated' in d.keys() for d in svs_config):
            return None
        # construct dict of form {chrom: [(start, end)]} Alu interval tuples -- will look for Alus in the self.overlap_events_dict
        alu_intervals = defaultdict(list)
        for elt_type, elt_list in self.overlap_events_dict.items():
            if 'Alu' in elt_type:
                for chrom, start, end in elt_list:
                    alu_intervals[chrom].append((int(start), int(end)))
        for chrom in alu_intervals.keys():
            random.shuffle(alu_intervals[chrom])
        # -------
        alu_pairs_dict = {}
        # construct the lists of viable flanking pairs for each SV config with a nonzero num_alu_mediated specified
        for sv_config in svs_config:
            sv_config_id = get_sv_config_identifier(sv_config)
            if 'num_alu_mediated' in sv_config.keys():
                # want a copy of the alu_intervals specific to this SV config so we can edit without changing the state of alu_intervals
                alu_intervals_copy = copy.deepcopy(alu_intervals)
                alu_pairs = []
                # recall: only DUPs and DELs may be made to be alu-mediated so length bounds will always be given as ints
                sv_min, sv_max = sv_config['min_length'], sv_config['max_length']
                # collecting twice as many Alu pairs as necessary to allow for collisions with previously-places SVs
                # --> terminate if we run out of viable matches (we're picking the left_alu without replacement, so stop when there
                # --> are no choices left)
                for _ in range(self.svtype_alu_mediated_counts[sv_config_id]):  # * 2): <-- just going to generate the exact number requested, see if it becomes an issue with collisions with previously-placed SVs
                    viable_matches = []
                    while len(viable_matches) == 0:
                        if len(alu_intervals_copy.keys()) == 0 or not any(len(v) > 1 for v in alu_intervals_copy.values()):
                            # break if there are no alus remaining to choose from
                            break
                        # choose a random starting (left) Alu (random across chromosome and position)
                        # --> only need to choose from chroms with at least 2 alus
                        rand_chrom = random.choice([chrom for chrom, alus in alu_intervals_copy.items() if len(alus) > 1])
                        alu1 = alu_intervals_copy[rand_chrom].pop()
                        # collect all the alus on the same chromosome that are within an appropriate dist (and repeat
                        # initial alu selection if there are no appropriate matches)
                        viable_matches = [ivl for ivl in alu_intervals_copy[rand_chrom] if sv_min <= np.abs(self.midpoint(*ivl) - self.midpoint(*alu1)) <= sv_max]
                    if len(viable_matches) > 0:
                        # need to pick this without replacement as well, right?
                        alu2 = random.choice(viable_matches)
                        alu_intervals_copy[rand_chrom].remove(alu2)
                        # -- update alu_intervals --
                        if len(alu_intervals_copy[rand_chrom]) == 0:
                            del alu_intervals_copy[rand_chrom]
                        # --------------------------
                        # need to also record chrom so we can yield that along with the interceding interval
                        alu_pairs.append((rand_chrom, alu1, alu2) if alu1[0] < alu2[0] else (rand_chrom, alu2, alu1))
                        # and only want to delete from the true alu_intervals dict the Alus that are actually used
                        # (not just the ones that didn't yield matches, because those might have matches under a different SV config's params)
                        # --> also need to remove from the general overlap_events_dict because these now cant be used for full overlap placement either
                        for alu in [alu1, alu2]:
                            alu_intervals[rand_chrom].remove(alu)
                            self.remove_alu_from_overlap_dict(rand_chrom, *alu)
                if len(alu_pairs) > 0:
                    random.shuffle(alu_pairs)
                    alu_pairs_dict[sv_config_id] = alu_pairs
        return alu_pairs_dict

    def get_alu_mediated_interval(self, sv_config_id):
        # # --- debug ---
        # print(f'get_alu_med_interval called\nBEFORE --> svtype_alu_mediated_counts = {self.svtype_alu_mediated_counts}')
        # # -------------
        # 0) if there are no remaining alu pairs to choose for this sv config, return None and zero out the remaining
        # counts in the svtype_alu...counts dict (so we don't keep looking in subsequent iterations)
        if len(self.alu_pairs[sv_config_id]) == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
            return None, None
        # 1) draw alu-mediated interval for a given DEL/DUP of given size range (don't need to give the size bounds
        # to this function because they were already used to assemble the list of viable alu pairs for this sv config)
        ivl_chrom, left_alu, right_alu = self.alu_pairs[sv_config_id].pop()
        # 2) construct the interceding interval (extending halfway through both Alus)
        ivl_start, ivl_end = (left_alu[1] + left_alu[0]) // 2, (right_alu[1] + right_alu[0]) // 2
        # 3) decrement the svtype_alu_mediated_counts entry for this SV config
        self.svtype_alu_mediated_counts[sv_config_id] -= 1
        # 3a) ... removing the dict entry if decremented to 0
        if self.svtype_alu_mediated_counts[sv_config_id] == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
        # # --- debug ---
        # print(f'get_alu_med_interval called\nAFTER --> svtype_alu_mediated_counts = {self.svtype_alu_mediated_counts}')
        # # -------------
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

    def __getitem__(self, sv_config_id, minsize, maxsize, elt_type=None, partial_overlap=False):
        # --> need to store input_elt_type to decrement the right entries at the end (could instead decrement
        # the count dict here but seems clearer to update both dicts together at the end)
        input_elt_type = 'ALL' if elt_type is None else elt_type
        if elt_type in self.overlap_events_dict.keys():  # <- elt_type given, and elements of matching type are in the dict
            # elt_type-specific branch: draw random element (recall: list is shuffled) of appropriate size
            rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None)
        elif len(self.overlap_events_dict.keys()) > 0:  # <- element dict is non-empty; elt_type not given, or is given as a prefix
            if elt_type is None:
                # to remove the chosen element from the right list in the dict we need a mapping from elt to elt_type
                elt_type_mapping = {elt: elt_type for (elt_type, elt_list) in self.overlap_events_dict.items() for elt in elt_list}
                # elt_type-agnostic branch: draw random element from combined list of elements across all types
                # --> ... and if allow_types is specified, then only want to grab elements in agreement with those
                rand_elt = next((elt for elt in itertools.chain.from_iterable(self.overlap_events_dict.values()) if
                                 (minsize <= self.get_intrvl_len(*elt) <= maxsize) and
                                 (self.elt_type_is_allowed(elt_type_mapping[elt]))), None)
            else:
                # elt_type prefix branch: choose an element whose name begins with elt_type
                # prefix_matches = [element_type for element_type in self.overlap_events_dict.keys() if element_type.startswith(elt_type)]
                prefix_matches = [element_type for element_type in self.overlap_events_dict.keys() if element_type.startswith(elt_type)]
                if len(prefix_matches) > 0:
                    random.shuffle(prefix_matches)
                    elt_type = prefix_matches.pop()
                    # make sure we choose an elt_type with a nonzero number of elements of valid size
                    while len(prefix_matches) > 0 and len([elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize]) == 0:
                        elt_type = prefix_matches.pop()
                    rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None)
                else:
                    rand_elt = None
        else:  # <- elt_type is NOT given and there AREN'T elements in the dict
            # don't assign an element (this SV will not overlap something)
            rand_elt = None
        if rand_elt is not None:
            # if we drew a non-None element, want to remove it from its list and decrement the count dictionary
            if elt_type is None:
                elt_type = elt_type_mapping[rand_elt]
            # TODO?: need to also make sure the event is removed from the available Alu intervals so we don't choose the same one and cause a collision
            #  --> not necessary if we always prioritize Alu-pair selection; is that necessarily bad? Seems potentially useful
            self.overlap_events_dict[elt_type].remove(rand_elt)
            if len(self.overlap_events_dict[elt_type]) == 0:
                del self.overlap_events_dict[elt_type]
            # -------------
            # decrement the relevant counts dict
            counts_dict = self.svtype_partial_overlap_counts if partial_overlap else self.svtype_overlap_counts
            counts_dict[sv_config_id][input_elt_type] -= 1
            # ... and when the number reaches 0, remove from dict so line 253 will not be triggered
            if counts_dict[sv_config_id][input_elt_type] == 0:
                del counts_dict[sv_config_id][input_elt_type]
                # ... and if the total counts for that svtype have been removed, remove the sv type dict entry
                if len(counts_dict[sv_config_id].keys()) == 0:
                    del counts_dict[sv_config_id]
        # check if we want to return a partially overlapping interval; if yes, compute one based on the elt chosen
        if partial_overlap:
            rand_elt = self.get_partially_overlapping_interval(*rand_elt, minsize, maxsize)
        # returning elt_type as well in order to have the retrieved type even when allow_types == 'ALL'
        return rand_elt, elt_type


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Special use case for utils.py: Commandline access to get_original_base_pos()')
    parser.add_argument('--genome_vcf', type=str, help='Path to vcf describing simulated SVs populating the synthetic genome')
    parser.add_argument('--query_loc', type=int, help='Query genome locus (given with respect to synthetic reference)')
    parser.add_argument('--query_chrom', type=int, help='Corresponding chromosome of the query locus')
    args = parser.parse_args()

    print('Input reference position corresponding to mutated reference position (%s) %d: %d' %
          (args.query_chrom, args.query_loc, get_original_base_pos(args.genome_vcf, args.query_loc, args.query_chrom)))

