from os import error
import random
import numpy as np
import argparse
from pysam import FastaFile
from processing import Error, FormatterIO, collect_args
from constants import Constants, Variant_Type, Symbols
#from arguments import collect_args
#import tracemalloc   # only for testing
import sys
from enum import Enum, unique
import time

# issues: insertion
# lowercase = invert
# _ = space
# ' = complement


time_start = time.time()

class ErrorDetection():
    def __init__(self):
        pass

    def is_overlapping(self, event_ranges, addition):
        # addition: tuple (start, end)
        # event_ranges: list containing tuples
        # checks if addition overlaps with any of the events already stored

        for event in event_ranges:
            if event[1] > addition[0] and event[0] < addition[1]:
                return True
        return False

    def any_overlapping(self, arr):
        for x, ele in enumerate(arr):
            if self.is_overlapping(arr[:x], ele):
                raise Exception("Overlapping Detected: {}".format(arr))
        
    def is_unique(self,transform):
        # source transformation must have unique symbols except dispersion events
        present = dict()
        for symbol in transform:
            if symbol in present:
                raise Exception("Source transformation {} does not have unique symbols!".format(transform))
            elif symbol != Symbols.DIS:
                present[symbol] = True

error_detection = ErrorDetection()

class Structural_Variant():
    def __init__(self, sv_type, length_ranges):
        '''
        sv_type: integer or tuple containing transformation (source, target)
        length_ranges: list containing tuple(s) (min_length, max_length)
        '''
        
        if isinstance(sv_type, int):
            self.type = Variant_Type(sv_type)
            self.source, self.target = Constants.SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.type = None
            self.source, self.target = sv_type
            self.name = str("".join(self.source)) + ">" + str("".join(self.target))
        
        error_detection.is_unique(self.source)
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
        return "<SV transformation {} -> {} taking up {} space>".format(self.source, self.target, self.req_space)
    
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

        decode_dict = {"invert": lambda string: complement(string[::-1]),
                       "identity": lambda string: string,
                       "complement": complement} 
        encoding = self.events_dict    # maps symbol like A or B to base pairs on reference 
        #print("Encode_dict: ", encoding)

        # find all blocks between dispersion events
        self.generate_blocks()
        
        changed_fragments = []
        for idx, block in enumerate(self.target_event_blocks):
            new_frag = ""
            for x, ele in enumerate(block):
                upper_str = ele.upper()         # changes all lowercase symbols (if there are any) to uppercase so we can map to the nucleotide bases 
    
                if any(c.islower() for c in ele):   # checks if lowercase symbols exist in ele
                    curr_piece = decode_dict["invert"](encoding[upper_str[0]].source_frag)   # remember decode_dict stores functions as the values
    
                elif upper_str[0] in encoding:
                    curr_piece = decode_dict["identity"](encoding[upper_str[0]].source_frag)
                
                elif ele.startswith(Symbols.DIS):               # _ refers to a space
                    raise Exception("Dispersion event detected within block: {}".format(self.target_event_blocks))
                else: 
                    raise Exception("Unknown {} symbol detected in target transformation {}".format(ele, self.target))
                
                new_frag += curr_piece
            first_event = encoding[self.source_event_blocks[idx][0].upper()]
            last_event = encoding[self.source_event_blocks[idx][-1].upper()]
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
        self.source_chr = ""
        self.source_frag = ""
        self.start = -1
        self.end = -1
    
    def __repr__(self):
        return "<Event {}>".format({"length": self.length, "symbol": self.symbol, "start": self.start, "end": self.end, "source_frag": self.source_frag})

class SV_Simulator():
    def __init__(self, ref_file, par_file):
        '''
        ref: fasta filename to reference genome
        svs: list (type, number, range of length of variant)

        '''

        global time_start
        print("Setting Up Simulator...")
        self.ref_file = ref_file
        self.ref_fasta = FastaFile(ref_file) 
        self.order_ids = self.ref_fasta.references
        self.len_dict = dict()
        for id in self.order_ids:
            self.len_dict[id] = self.ref_fasta.get_reference_length(id)
        print("Total base count: ", sum(self.ref_fasta.lengths))

        self.formatter = FormatterIO(ref_file)
        svs_config = self.formatter.yaml_to_var_list(par_file)
        self.initialize_svs(svs_config)

        print("Finished Setting up Simulator in {} seconds\n".format(time.time() - time_start))
        time_start = time.time()
    
    def __repr__(self):
        message = "SVs1: " + str([self.svs[x].name for x in range(len(self.svs)) if self.svs[x].hap[0]]) + "\n"
        message += "SVs2: " + str([self.svs[y].name for y in range(len(self.svs)) if self.svs[y].hap[1]]) + "\n"
        return message
    
    def initialize_svs(self, svs_config):
        self.svs = []
        for sv_config in svs_config:
            for num in range(sv_config[1]):
                sv = Structural_Variant(sv_config[0], sv_config[2]) # inputs: SV type, range of lengths
                draw = random.randint(1,4)
                if draw == 3 or draw == 4:   # sv applies to both haplotypes
                    sv.ishomozygous = 1
                    sv.hap = [True, True]
                elif draw == 2:
                    sv.ishomozygous = 0
                    sv.hap = [False, True]
                elif draw == 1:
                    sv.ishomozygous = 0
                    sv.hap = [True, False]
                self.svs.append(sv) 
    
    def reinitialize_svs(self):
        for sv in self.svs:
            sv.active = False

    def produce_variant_genome(self, fasta1_out, fasta2_out, bedfile, initial_reset = True, verbose = False):
        '''
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        '''
        global time_start
        if initial_reset:
            self.formatter.reset_file(fasta1_out)
            self.formatter.reset_file(fasta2_out)

        # edit chromosome
        ref_fasta = self.ref_fasta
        self.rand_edit_svs(ref_fasta)

        # organize edits and export
        active_svs = [sv for sv in self.svs if sv.active]

        for x in range(2):
            edits_dict = dict()
            for id in self.order_ids:
                edits_dict[id] = []

            if x == 0:
                fasta_out = fasta1_out
            elif x == 1:
                fasta_out = fasta2_out
            for sv in active_svs:
                if sv.hap[x]:
                    for frag in sv.changed_fragments:
                        edits_dict[frag[0]].append(frag[1:])

            for id in self.order_ids:
                # account for homozygous and heterogeneous variants
                edits_x = edits_dict[id]
                #print("Edits_x: ", edits_x)
                error_detection.any_overlapping(edits_x)
                
                # export edited chromosomes to FASTA files
                self.formatter.export_piece(id, edits_x, fasta_out, ref_fasta, verbose = verbose)

                print("ID {} altered and saved in fasta file {} in {} seconds\n".format(id, fasta_out, time.time() - time_start))
                time_start = time.time()

        # export variant data to BED file
        self.formatter.export_to_bedpe(active_svs, id, bedfile, reset_file = initial_reset)

        initial_reset = False

        return True
    
    def rand_select_svs(self, svs, ref_fasta):
        # randomly position SVs and store reference fragments

        def percent_N(seq):
            total = 0
            if len(seq) == 0:  # avoid ZeroDivisionError
                return 0
            for char in seq:
                if char == "N":
                    total += 1
            return total / len(seq)
        def is_overlapping(event_ranges, addition):
            # addition: tuple (start, end)
            # event_ranges: list containing tuples
            # checks if addition overlaps with any of the events already stored

            for event in event_ranges:
                if event[1] > addition[0] and event[0] < addition[1]:
                    return True
            return False
        def generate_seq(length):
            # helper function for insertions
            # generates random sequence of bases of given length
            rand_seq = ""
            base_map = {1:"A", 2: "T", 3: "G", 4: "C"}
            for x in range(length):
                rand_seq += base_map[random.randint(1,4)]
            return rand_seq
        def get_rand_chr():
            # random assignment of SV to a chromosome
            rand_id = self.order_ids[random.randint(0, len(self.order_ids)-1)]
            chr_len = self.len_dict[rand_id]

            return rand_id, chr_len
        
        self.reinitialize_svs()
        
        self.event_ranges = []
        for sv in svs:
            #print("Current SV: ", sv)
            tries = 0 # number of attempts to place sv randomly
            valid = False
            while not valid:
                tries += 1
                valid = True
                
                if tries > 100:
                    print("Failure to simulate SV \"{}\"".format(sv))
                    valid = False
                    break
                
                rand_id, chr_len = get_rand_chr()
                if chr_len - sv.req_space - 1 <= 0:
                    raise Exception("{} size is too big for chromosome!".format(sv))
                else:
                    start_pos = random.randint(0, chr_len - sv.req_space - 1)
                    for sv_event in sv.source_events:
                        #print("Symbol: {}, Valid: {}".format(sv_event.symbol, valid))

                        # store start and end position and reference fragment
                        sv_event.start, sv_event.end = start_pos, start_pos + sv_event.length
                        sv_event.source_chr = rand_id
                        frag = ref_fasta.fetch(rand_id, sv_event.start, sv_event.end)
                        sv_event.source_frag = frag
                        start_pos += sv_event.length

                        if sv_event.symbol.startswith(Symbols.DIS):
                            continue

                        # check to see if chosen spot is a valid position
                        #print("Percent_N, {} -> {}".format(frag, percent_N(frag)))
                        #print("Is_overlapping, {} -> {}".format((self.event_ranges, (sv_event.start, sv_event.end)), is_overlapping(self.event_ranges, (sv_event.start, sv_event.end))))
                        if percent_N(frag) > 0.05 or is_overlapping(self.event_ranges, (sv_event.start, sv_event.end)):
                            valid = False

            
            # dispersion events (indicated with a space) are not off-limits to more events
            if valid:
                sv.active = True
                self.event_ranges.extend([(event.start, event.end) for event in sv.source_events if not event.symbol.startswith(Symbols.DIS)])
                sv.start = sv.source_events[0].start
                sv.end = sv.source_events[-1].end

                # handles insertions - these event symbols only show up in target transformation
                for event in sv.events_dict.values():
                    if event.source_frag == "" and event.length > 0:
                        event.source_frag = generate_seq(event.length)
                
                # error detection
                error_detection.any_overlapping(self.event_ranges)
                
                #print("Events Dict after random positions selected: ", sv.events_dict)
                #print("\n")
                
        return self.event_ranges

    def rand_edit_svs(self, ref_fasta):
        # select random positions for SVs
        self.rand_select_svs(self.svs, ref_fasta)

        for sv in self.svs:
            if sv.active:
                # make edits and store in sv object
                sv.change_fragment()
        

if __name__ == "__main__":

    args = collect_args()
    #tracemalloc.start()
    fasta_in = args[0]
    yaml_in = args[1]
    fasta1_out = args[2]
    fasta2_out = args[3]
    bed_out = args[4]

    '''fasta_in = "debugging/inputs/test.fna"
    yaml_in = "debugging/inputs/par_test.yaml"
    #test_svs = [[12,[(30,100)]], [16,[(300,8000)]],[8,[(500,1000)]]]
    fasta1_out = "debugging/inputs/test1_out.fna"
    fasta2_out = "debugging/inputs/test2_out.fna"
    bed_out = "debugging/inputs/out.bed"'''

    sim = SV_Simulator(fasta_in, yaml_in)
    print(sim.svs)
    #print(sim.rand_edit_svs("Chromosome19", sim.ref_fasta))
    sim.produce_variant_genome(fasta1_out, fasta2_out, bed_out, verbose = False)
    print("\n" + str(sim))

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    #tracemalloc.stop()


