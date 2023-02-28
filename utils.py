from constants import *
import os
import random

def is_overlapping(event_ranges, addition):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored

    # *** these need to be strict inequalities to allow for abutting insertions/deletions
    for event in event_ranges:
        # *** ... and with strict inequalities these are all given by the same expression
        # if addition[1] == addition[0]: # addition is an insertion
        #     if addition[0] > event[0] and addition[1] < event[1]:
        #         print('is_overlapping: case 1')
        #         return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])
        # if event[0] == event[1]:  # event is insertion
        #     if event[0] > addition[0] and event[0] < addition[1]:
        #         print('is_overlapping: case 2')
        #         print(f'event={event}; addition={addition}')
        #         return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])
        if event[1] > addition[0] and event[0] < addition[1]:
            print('is_overlapping: case 3')
            return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])

    return False

def fail_if_any_overlapping(arr):
    # will raise Exception if any overlap between intervals is found
    # arr: list of tuples
    # debug
    print(f'overlap check on {arr}')
    for x, ele in enumerate(arr):
        print(f'{x}, {ele}')
        if is_overlapping(arr[:x], ele):
            raise Exception(is_overlapping(arr[:x], ele)[1])

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

        if Symbols.DUP_MARKING.value in symbol:
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

def parse_bed_file(bed_fname):
    """
    reads bed file (intended use: processing repeatmasker elements to be considered for randomized
    event overlap) and returns list of (chr, start, end) tuples representing the intervals of each event in the file
    --> logic taken from bed_iter() and parse_bed_line() in cue (seq/io.py)
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
            intervals_list.append((chr_name, start, end))
    return intervals_list
