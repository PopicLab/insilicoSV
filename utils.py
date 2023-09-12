import os
import random
import itertools
from constants import *
from collections import defaultdict


# NestedDict helper class
class NestedDict(defaultdict):
    """Nested defaultdict"""

    def __call__(self):
        return NestedDict(self.default_factory)


def is_overlapping(event_ranges, addition, called_from_helper=False):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored
    for event in event_ranges:
        if event[1] > addition[0] and event[0] < addition[1]:
            if called_from_helper:
                raise Exception("Overlap between {} and {}".format(event[0:2], addition[0:2]))
            else:
                return True

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
    base_map = {1:"A", 2: "T", 3: "G", 4: "C"}
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
    max_str = str(sv_config['min_length']) if isinstance(sv_config['min_length'], int) else '_'.join(map(str, sv_config['min_length']))
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
        self.svtype_alu_mediated_counts = defaultdict(NestedDict(int))
        self.get_num_overlap_counts(config)

        # overlap_events can either be given as a single .bed file or multiple
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
            if 'num_overlap' in sv.keys():
                if isinstance(sv['num_overlap'], int):
                    # if num_overlap given as a singleton integer, label type as 'ALL' if zero or multiple allow_types
                    # are given, otherwise set type to the single allow_type provided
                    if self.allow_types is not None and len(self.allow_types) == 1:
                        self.svtype_overlap_counts[sv_config_key][self.allow_types[0]] = sv['num_overlap']
                    else:
                        self.svtype_overlap_counts[sv_config_key]['ALL'] = sv['num_overlap']
                else:
                    # set up correspondence with the element types as they were given in the config file
                    for num_ovlp in sv['num_overlap']:
                        self.svtype_overlap_counts[sv_config_key][self.allow_types[sv['num_overlap'].index(num_ovlp)]] = num_ovlp
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
                if (allow_chroms and chr_name not in allow_chroms) or (allow_types and not any([elt_name in elt_type for elt_name in allow_types])):
                    continue
                self.overlap_events_dict[elt_type].append((chr_name, start, end))
        # shuffle the lists of each element type for random selection
        for k in self.overlap_events_dict.keys():
            random.shuffle(self.overlap_events_dict[k])

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
        # -------
        alu_pairs_dict = {}
        # construct the lists of viable flanking pairs for each SV config with a nonzero num_alu_mediated specified
        for sv_config in svs_config:
            sv_config_id = get_sv_config_identifier(sv_config)
            if 'num_alu_mediated' in sv_config.keys():
                alu_pairs = []
                # recall: only DUPs and DELs may be made to be alu-mediated so length bounds will always be given as ints
                sv_min, sv_max = sv_config['min_length'], sv_config['max_length']
                # collecting twice as many Alu pairs as necessary to allow for collisions with previously-places SVs
                for _ in range(self.svtype_alu_mediated_counts[sv_config_id] * 2):
                    viable_matches = []
                    while len(viable_matches) == 0:
                        # choose a random starting (left) Alu (random across chromosome and position)
                        rand_chrom = random.choice(list(alu_intervals.keys()))
                        left_alu = random.choice(alu_intervals[rand_chrom])
                        # collect all the alus on the same chromosome that are within an appropriate dist (and repeat
                        # initial alu selection if there are no appropriate matches)
                        viable_matches = [ivl for ivl in alu_intervals[rand_chrom] if sv_min <= ivl[0] - left_alu[1] <= sv_max]
                    right_alu = random.choice(viable_matches)
                    # need to also record chrom so we can yield that along with the interceding interval
                    alu_pairs.append((rand_chrom, left_alu, right_alu))
                random.shuffle(alu_pairs)
                alu_pairs_dict[sv_config_id] = alu_pairs
        return alu_pairs_dict

    def get_alu_mediated_interval(self, sv_config_id):
        # 1) draw alu-mediated interval for a given DEL/DUP of given size range (don't need to give the size bounds
        # to this function because they were already used to assemble the list of viable alu pairs for this sv config)
        # TODO: add a test case for addressing collisions between extracted Alu-mediated intervals and already-placed SVs
        ivl_chrom, left_alu, right_alu = random.choice(self.alu_pairs[sv_config_id])
        # 2) construct the interceding interval (extending halfway through both Alus)
        ivl_start, ivl_end = (left_alu[1] + left_alu[0]) // 2, (right_alu[1] + right_alu[0]) // 2
        # 3) decrement the svtype_alu_mediated_counts entry for this SV config
        self.svtype_alu_mediated_counts[sv_config_id] -= 1
        # 3a) ... removing the dict entry if decremented to 0
        if self.svtype_alu_mediated_counts[sv_config_id] == 0:
            del self.svtype_alu_mediated_counts[sv_config_id]
        return ivl_chrom, ivl_start, ivl_end, 'ALU_MEDIATED'  # <- ALU_MEDIATED to be used as retrieved_type

    @staticmethod
    def get_intrvl_len(chr, st, end):
        return int(end) - int(st)

    def __getitem__(self, sv_config_id, minsize, maxsize, elt_type=None):
        # # debug
        # print(f'getitem() called with elt_type = {elt_type}')
        # print(f'svtype_overlap_counts = {self.svtype_overlap_counts}')
        # --> need to store input_elt_type to decrement the right entries at the end (could instead decrement
        # the count dict here but seems cleared to udpate both dicts together at the end)
        input_elt_type = elt_type
        if elt_type in self.overlap_events_dict.keys():  # <- elt_type given, and elements of matching type are in the dict
            # elt_type-specific branch: draw random element (recall: list is shuffled) of appropriate size
            rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None)
        elif len(self.overlap_events_dict.keys()) > 0:  # <- element dict is non-empty; elt_type not given, or is given as a prefix
            if elt_type is None:
                # elt_type-agnostic branch: draw random element from combined list of elements across all types
                rand_elt = next((elt for elt in itertools.chain.from_iterable(self.overlap_events_dict.values()) if minsize <= self.get_intrvl_len(*elt) <= maxsize), None)
                # to remove the chosen element from the right list in the dict we need a mapping from elt to elt_type
                elt_type_mapping = {elt: elt_type for (elt_type, elt_list) in self.overlap_events_dict.items() for elt in elt_list}
            else:
                # elt_type prefix branch: choose an element whose name begins with elt_type
                prefix_matches = [element_type for element_type in self.overlap_events_dict.keys() if element_type.startswith(elt_type)]
                if len(prefix_matches) > 0:
                    elt_type = random.sample(prefix_matches, 1)[0]
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
            self.overlap_events_dict[elt_type].remove(rand_elt)
            if len(self.overlap_events_dict[elt_type]) == 0:
                del self.overlap_events_dict[elt_type]
            # -------------
            # decrement the self.svtype_overlap_counts dict
            self.svtype_overlap_counts[sv_config_id][input_elt_type] -= 1
            # ... and when the number reaches 0, remove from dict so line 253 will not be triggered
            if self.svtype_overlap_counts[sv_config_id][input_elt_type] == 0:
                del self.svtype_overlap_counts[sv_config_id][input_elt_type]
                # ... and if the total counts for that svtype have been removed, remove the sv type dict entry
                if len(self.svtype_overlap_counts[sv_config_id].keys()) == 0:
                    del self.svtype_overlap_counts[sv_config_id]
        # returning elt_type as well in order to have the retrieved type even when allow_types == 'ALL'
        return rand_elt, elt_type

