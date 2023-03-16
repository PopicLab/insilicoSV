from constants import *
import os
import random

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
        # instead of this, moving the exception to be raised by is_overlapping()
        is_overlapping(arr[:x], ele, called_from_helper=True)
        # if is_overlapping(arr[:x], ele, called_from_helper=True):
        #     raise Exception(is_overlapping(arr[:x], ele)[1])

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
        elif symbol != Symbols.DIS_MARKING.value:   # exclude dispersion events because they will always appear the same for user inputs
            present[symbol] = True

        if Symbols.DUP_MARKING.value in symbol:
            raise Exception("Duplication marking (') not allowed in source sequence {}".format(source))
        if any(c.islower() for c in symbol):
            raise Exception("Only original symbols may appear in source sequence, and they must also be uppercase: {}".format(source))
    
    assert sum([1 for symbol in source if symbol.startswith(Symbols.DIS_MARKING.value)]) == sum([1 for symbol in target if symbol.startswith(Symbols.DIS_MARKING.value)]), "Number of dispersion events must be equal between source and target sequence"

def remove_file(file):
    # removes file if it exists
    if os.path.exists(file):
        os.remove(file)

def reset_file(filename):
    #print("Overwritting File {}...".format(filename))
    with open(filename, "w") as f_reset:
        f_reset.truncate()

def generate_seq(length, random_gen=None):
    # helper function for insertions
    # generates random sequence of bases of given length
    base_map = {1:"A", 2: "T", 3: "G", 4: "C"}
    return ''.join([base_map[random_gen.randint(1, 4)] for _ in range(length)])

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
    # function to create slightly mutated version of an input sequence (to be used to create
    # the quasi-repeated section of a divergent repeat/divergent dDUP)
    # --> p given as the probability of changing the base, chosen from U(0.5,1.0)
    p = random.uniform(0.5,1.0)
    return ''.join([b if random.random() > p else random.choice(list({"A", "C", "T", "G"} - {b})) for b in seq.upper()])

def parse_bed_file(bed_fname, keep_type=False):
    """
    reads bed file (intended use: processing repeatmasker elements to be considered for randomized
    event overlap) and returns list of (chr, start, end) tuples representing the intervals of each event in the file
    --> logic taken from bed_iter() and parse_bed_line() in cue (seq/io.py)
    - keep_type: optional flag to extract intervals with the fourth column string that in the case of a repeatmasker
                    bed file will give the repetitive element type
    """
    intervals_list = []
    with open(bed_fname, 'r') as bed_file:
        for line in bed_file:
            if line.startswith('#') or line.isspace():
                continue
            # parse line
            fields = line.strip().split()
            assert len(fields) >= 3, "Unexpected number of fields in BED: %s" % line
            chr_name, start, end = fields[:3]
            if keep_type:
                intervals_list.append((chr_name, start, end, fields[3]))
            else:
                intervals_list.append((chr_name, start, end))
    return intervals_list

def process_overlap_events(config):
    overlap_events = []
    if type(config.overlap_events['bed']) is list:
        overlap_events = []
        for bed_path in config.overlap_events['bed']:
            # need to extract bed record intervals with the element type given in column 4 (keep_type=True)
            overlap_events.extend(parse_bed_file(bed_path, keep_type=True))
    else:
        overlap_events = parse_bed_file(config.overlap_events['bed'], keep_type=True)
    random.shuffle(overlap_events)
    # filter on allowed repetitive element types (if any are given)
    # --> if given, allow_types must be given as a list of strings
    if 'allow_types' in config.overlap_events.keys():
        # remove the repetitive element type at the filtering step (we just care that it's an allowed type)
        overlap_events = [ev[:-1] for ev in overlap_events if ev[-1] in config.overlap_events['allow_types']]
    return overlap_events
