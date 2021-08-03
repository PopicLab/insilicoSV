
from constants import *
from processing import ErrorDetection
import random

class Structural_Variant():
    def __init__(self, sv_type, length_ranges, source = None, target = None):
        '''
        Initializes SV's transformation and sets up its events, along with several other basic attributes like zygosity

        sv_type: Enum either specifying one of the prewritten classes or a Custom transformation, in which case source and target are required
        length_ranges: list containing tuple(s) (min_length, max_length)
        source: tuple representing source sequence, optional
        target: tuple representing target sequence, optional
        '''

        self.type = sv_type
        if self.type != Variant_Type.Custom:
            self.source, self.target = SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.source, self.target = source, target
            self.name = str("".join(self.source)) + ">" + str("".join(self.target))
        
        # events like dispersions will appear as the same symbol, so it's important to add unique tags to differentiate them
        ErrorDetection.validate_symbols(self.source)
        self.source_unique_char, self.target_unique_char = self.add_unique_ids(self.source), self.add_unique_ids(self.target)

        # initialize event classes
        self.req_space = -1
        self.source_events = []
        self.events_dict = dict()
        self.initialize_events(length_ranges)
        #print("Events Dict: ", self.events_dict)
        self.source_event_blocks = []
        self.target_event_blocks = []
        self.target_block_chrs = []   # stores the assigned target chr for each block to export to BEDPE file

        # specifies if sv is unable to be simulated due to random placement issues
        # will be turned on later
        self.active = False

        # 1 = homozygous, 0 = heterozygous
        self.ishomozygous = Zygosity.UNDEFINED

        # stores list of booleans specifying if SV will be applied to certain haplotype for assigned chromosome
        self.hap = [False, False]
    
    def __repr__(self):
        return "<SV transformation {} -> {} taking up {} non-dispersion spaces>".format(''.join(self.source), ''.join(self.target), sum([event.length for event in self.source_events if not event.symbol.startswith(Symbols.DIS)]))
    
    def add_unique_ids(self, transformation):
        # if dispersion events exist in transformation, tag on unique ids to make them distinct as they all are "_"
        # unique ids later necessary to map symbol to event
        unique_transform = []
        unique_id = 1
        for component in transformation:
            if component != Symbols.DIS:
                unique_transform.append(component)
            else:
                unique_transform.append(component + str(unique_id))
                unique_id += 1
        return tuple(unique_transform)

    def initialize_events(self, lengths):
        '''
        Initializes event classes and creates a mapping of symbol to event

        lengths: list of tuples specifying min and max length for events within SV
        -> returns list of events in source sequence
        '''
        # collect all unique symbols present in both source and target sequences - include target as there may be insertions
        # note that the symbols represent events
        all_symbols = []
        for ele in self.source_unique_char + self.target_unique_char:
            if ele.upper() not in all_symbols:
                all_symbols.append(ele.upper())
        all_symbols.sort()  # user inputs symbol lengths in lexicographical order

        # symbols_dict: (key = symbol, value = (chosen length, length range))
        # determine length of events/symbols
        symbols_dict = dict()
        if len(lengths) > 1:    # values given by user represents custom ranges for each event symbol of variant in lexicographical order 
            assert (len(lengths) == len(all_symbols)) 
            for idx, symbol in enumerate(all_symbols):
                symbols_dict[symbol] = (random.randint(lengths[idx][0], lengths[idx][1]), lengths[idx])

        elif len(lengths) == 1: # value given by user represents length (same range) of each event within variant in lexicographical order
            for symbol in all_symbols:
                symbols_dict[symbol] = (random.randint(lengths[0][0], lengths[0][1]), lengths[0])

        else:
            raise Exception("Lengths parameter expects at least one tuple")
        symbols_dict[Symbols.PLACEHOLDER] = (0, (0,0))

        # initialize event classes
        for idx, symbol in enumerate(all_symbols):
            event = Event(self, symbols_dict[symbol][0], symbols_dict[symbol][1], symbol)
            self.events_dict[symbol] = event
        
        for symbol in self.source_unique_char:
            self.source_events.append(self.events_dict[symbol])
        
        self.req_space = sum([event.length for event in self.source_events])
        
        return self.source_events
    
    def generate_blocks(self):
        '''
        Groups together source and target symbols between dispersion events (_)
        Each block of symbols belongs to the same target chromosome
        change_fragment outputs edits using block's start and end position

        -> returns list of lists
        '''
        def find_blocks(transformation):
            # transformation: tuple of strings
            # -> returns list of lists of strings
            # Ex. ("A","B","_","C","D") -> [["A","B"], ["C","D"]]
            # Ex. ("A","B","_","_") -> [["A","B"],[],[]]
            blocks = [[]]
            for symbol in transformation:
                if not symbol.startswith(Symbols.DIS):
                    blocks[-1].append(symbol)
                else:
                    blocks.append([])
            return blocks
        self.source_event_blocks = find_blocks(self.source_unique_char)
        self.target_event_blocks = find_blocks(self.target_unique_char)

        return self.target_event_blocks

    def change_fragment(self):

        def complement(bases):
            # bases: str
            output = ""
            base_complements = {"A": "T", "T": "A", "G": "C", "C": "G", "N":"N"}
            for base in bases.upper():
                if base in base_complements:
                    output += base_complements[base]
                else:
                    output += base
                    print("Error: Unknown base \'{}\' detected".format(base))
            
            return output

        decode_funcs = {"invert": lambda string: complement(string[::-1]),
                       "identity": lambda string: string,
                       "complement": complement} 
        encoding = self.events_dict    # maps symbol like A or B to base pairs on reference 
        #print("Encode_dict: ", encoding)

        # find all blocks of symbols between dispersion events
        self.generate_blocks()

        changed_fragments = []
        for idx, block in enumerate(self.target_event_blocks):
            new_frag = ""
            for x, ele in enumerate(block):
                upper_str = ele.upper()         # used to find corresponding event from encoding, all keys in encoding are in uppercase
    
                if any(c.islower() for c in ele):   # checks if lowercase symbols exist in ele, represents an inversion
                    curr_piece = decode_funcs["invert"](encoding[upper_str].source_frag)   
    
                elif upper_str in encoding:   # take original fragment, no changes
                    curr_piece = decode_funcs["identity"](encoding[upper_str].source_frag)
                
                elif ele.startswith(Symbols.DIS):  # DIS = dispersion event ("_")
                    raise Exception("Dispersion event detected within block: {}".format(self.target_event_blocks))
                else: 
                    raise Exception("Unknown {} symbol detected in target transformation {}".format(ele, self.target))  # all possible symbols should have been logged and mapped to event
                
                new_frag += curr_piece
            
            # new_frag will replace events from the source sequence, so imagine from the point of view of the source sequence
            first_event = encoding[self.source_event_blocks[idx][0].upper()]
            last_event = encoding[self.source_event_blocks[idx][-1].upper()]

            # all events in source sequence MUST have a valid start & end position and a valid source chromosome
            assert (first_event.start != -1 and last_event.end != -1 and first_event.source_chr != None)
            changed_fragments.append([first_event.source_chr, first_event.start, last_event.end, new_frag])
            self.target_block_chrs.append(first_event.source_chr)  # used for exporting to bed file
                
        #print("ref {} -> transformed {} for transformation {}".format(ref_piece, change_genome, self.type.value))
        self.changed_fragments = changed_fragments
        print("Target Chromosomes: ", self.target_block_chrs)

        return changed_fragments

class Event():
    '''represents the symbols, also known as the "events," within a SV transformation'''
    def __init__(self, sv_parent, length, length_range, symbol):
        '''
        sv_parent: Structural Variant, event is always part of larger SV
        '''
        self.sv_parent = sv_parent
        self.length = length
        self.length_range = length_range
        self.symbol = symbol # refers to symbol in SV's transformation
        self.source_chr = None
        self.source_frag = ""
        self.start = -1
        self.end = -1
    
    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end, "source_frag": self.source_frag,
                        "source_chr": self.source_chr, "target_chr": self.target_chr})

