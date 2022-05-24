from constants import *
import os
import random

def is_overlapping(event_ranges, addition):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored

    for event in event_ranges:
        if addition[1] == addition[0]: # addition is an insertion
            if addition[0] >= event[0] and addition[1] <= event[1]: 
                return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])
        if event[0] == event[1]:  # event is insertion
            if event[0] >= addition[0] and event[0] <= addition[1]:
                return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])
        if event[1] > addition[0] and event[0] < addition[1]:
            return True, "Overlap between {} and {}".format(event[0:2], addition[0:2])

    return False

def fail_if_any_overlapping(arr):
    # will raise Exception if any overlap between intervals is found
    # arr: list of tuples
    for x, ele in enumerate(arr):
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
    if len(seq) == 0:  # avoid ZeroDivisionError
        return 0
    else:
        return seq.count('N') / len(seq)


def complement(seq):
    # seq: str
    output = ""
    base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    for base in seq.upper():
        if base in base_complements:
            output += base_complements[base]
        else:
            raise ValueError("Unknown base \'{}\' detected in reference, complement of base not taken".format(base))

    return output
