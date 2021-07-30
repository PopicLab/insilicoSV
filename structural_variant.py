
from constants import *
from processing import ErrorDetection
import random

class Structural_Variant():
    def __init__(self, sv_type, length_ranges):
        '''
        sv_type: string or tuple containing transformation (source, target)
        length_ranges: list containing tuple(s) (min_length, max_length)
        '''
        
        if isinstance(sv_type, str):
            self.type = Variant_Type(sv_type)
            self.source, self.target = SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.type = None
            self.source, self.target = sv_type
            self.name = str("".join(self.source)) + ">" + str("".join(self.target))
        
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

        # in the event that sv is unable to be simulated due to random placement issues, will be turned on later
        self.active = False

        # 1 = homozygous, 0 = heterozygous
        self.ishomozygous = -1
        self.hap1 = -1
        self.hap2 = -1
    
    def __repr__(self):
        return "<SV transformation {} -> {} taking up {} space>".format(''.join(self.source), ''.join(self.target), self.req_space)
    
    def add_unique_ids(self, transformation):
        # if dispersion events exist in transformation, tag on unique ids to make them distinct as they all are "_"
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

        # collect all unique symbols to account for insertions
        all_symbols = []
        for ele in self.source_unique_char + self.target_unique_char:
            if ele.upper() not in all_symbols:
                all_symbols.append(ele.upper())
        all_symbols.sort()

        # symbols_dict: (key = symbol, value = length)
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
        # groups together symbols between dispersion events (_)
        # returns list of lists
        # Ex. AB_CD -> [["A","B"], ["C","D"]]
        def find_blocks(transformation):
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
            output = ""
            for base in bases.upper():
                if base == "A":
                    output += "T"
                elif base == "T":
                    output += "A"
                elif base == "G":
                    output += "C"
                elif base == "C":
                    output += "G"
                elif base == "N":
                    output += "N"
                else:
                    output += base
                    print("Error: Unknown base \'{}\' detected".format(base))
            
            return output

        decode_funcs = {"invert": lambda string: complement(string[::-1]),
                       "identity": lambda string: string,
                       "complement": complement} 
        encoding = self.events_dict    # maps symbol like A or B to base pairs on reference 
        #print("Encode_dict: ", encoding)

        # find all blocks between dispersion events
        self.generate_blocks()
        
        changed_fragments = []
        for idx, block in enumerate(self.target_event_blocks):
            new_frag = ""
            # new_frag will replace events from the source sequence
            # imagine from the point of view of the source sequence
            first_event = encoding[self.source_event_blocks[idx][0].upper()]
            last_event = encoding[self.source_event_blocks[idx][-1].upper()]
            for x, ele in enumerate(block):
                upper_str = ele.upper()         # changes all lowercase symbols (if there are any) to uppercase so we can map to the nucleotide bases 
    
                if any(c.islower() for c in ele):   # checks if lowercase symbols exist in ele
                    curr_piece = decode_funcs["invert"](encoding[upper_str].source_frag)   
    
                elif upper_str in encoding:
                    curr_piece = decode_funcs["identity"](encoding[upper_str].source_frag)
                
                elif ele.startswith(Symbols.DIS):               # _ refers to a space
                    raise Exception("Dispersion event detected within block: {}".format(self.target_event_blocks))
                else: 
                    raise Exception("Unknown {} symbol detected in target transformation {}".format(ele, self.target))
                
                new_frag += curr_piece
                encoding[upper_str].target_chr = first_event.source_chr   # source symbols being replaced, so their source_chr is the target
            changed_fragments.append([first_event.source_chr, first_event.start, last_event.end, new_frag])
                
        #print("ref {} -> transformed {} for transformation {}".format(ref_piece, change_genome, self.type.value))
        self.changed_fragments = changed_fragments

        return changed_fragments

class Event():
    '''represents the symbols/events within a SV transformation'''
    def __init__(self, sv_parent, length, length_range, symbol):
        '''
        sv_parent: Structural Variant, event is a part of larger SV
        '''
        self.sv_parent = sv_parent
        self.length = length
        self.length_range = length_range
        self.symbol = symbol # refers to symbol in SV's transformation
        self.source_chr = "None"
        self.source_frag = ""
        self.target_chr = "None"
        self.start = -1
        self.end = -1
    
    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end, "source_frag": self.source_frag,
                        "source_chr": self.source_chr, "target_chr": self.target_chr})

