from constants import *
import os

def is_overlapping(event_ranges, addition):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored

    for event in event_ranges:
        if addition[0] == addition[1] and addition[1] == event[0] and event[0] == event[1]:  # reject insertions at the same position
            return True
        if event[1] > addition[0] and event[0] < addition[1]:
            return True

    return False

def fail_if_any_overlapping(arr):
    # will raise Exception if any overlap between intervals is found
    # arr: list of tuples
    for x, ele in enumerate(arr):
        if is_overlapping(arr[:x], ele):
            raise Exception("Overlapping Detected: {}".format(arr))

def validate_symbols(transform):
    '''
    Ensures that transform has unique symbols other than dispersion events - specifically used for source sequence

    transform: str, source sequence
    '''
    present = dict()
    for symbol in transform:
        if symbol in present:
            raise Exception("Source transformation {} does not have unique symbols!".format(transform))
        elif symbol != Symbols.DIS.value:   # exclude dispersion events because they will always appear the same for user inputs
            present[symbol] = True

        if Symbols.DUP_MARKING.value in symbol:
            raise Exception("Duplication marking (') not allowed in source sequence {}".format(transform))
        if any(c.islower() for c in symbol):
            raise Exception("Only original symbols may appear in source sequence, and they must also be uppercase: {}".format(transform))

def remove_file(file):
    # removes file if it exists
    if os.path.exists(file):
        os.remove(file)


