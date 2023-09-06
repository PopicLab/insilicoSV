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
        # dict containing the counts of overlap for each SV type and known element type
        self.svtype_overlap_counts = defaultdict(NestedDict(int))
        self.get_num_overlap_counts(config)

        # overlap_events can either be given as a single .bed file or multiple
        if type(config['overlap_events']['bed']) is list:
            for bed_path in config['overlap_events']['bed']:
                self.parse_bed_file(bed_path, allow_chroms=self.allow_chroms, allow_types=self.allow_types)
        else:
            self.parse_bed_file(config['overlap_events']['bed'], allow_chroms=self.allow_chroms, allow_types=self.allow_types)

    def get_num_overlap_counts(self, config):
        # populate nested dict of the form {sv type: {element type: num_overlap}}
        for sv in config['SVs']:
            if isinstance(sv['num_overlap'], int):
                self.svtype_overlap_counts[sv['type']]['ALL'] = sv['num_overlap']
            else:
                # set up correspondence with the element types as they were given in the config file
                for num_ovlp in sv['num_overlap']:
                    self.svtype_overlap_counts[sv['type']][self.allow_types[sv['num_overlap'].index(num_ovlp)]] = num_ovlp

    def parse_bed_file(self, bed_fname, allow_chroms=None, allow_types=None):
        """
        reads bed file (intended use: processing repeatmasker elements to be considered for randomized
        event overlap) and addes (chr, start, end) tuples representing the intervals to a dict keyed on elt type
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

    @staticmethod
    def get_intrvl_len(chr, st, end):
        return int(end) - int(st)

    def __getitem__(self, minsize, maxsize, elt_type=None):
        if elt_type in self.overlap_events_dict.keys():  # <- elt_type given, and elements of matching type are in the dict
            # elt_type-specific branch: draw random element (recall: list is shuffled) of appropriate size
            rand_elt = next((elt for elt in self.overlap_events_dict[elt_type] if minsize <= self.get_intrvl_len(*elt) <= maxsize), None)
        # elif elt_type is None and len(self.overlap_events_dict.keys()) > 0:  # <- elt_type not given, and the element dict is non-empty
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
        return rand_elt

