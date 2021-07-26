import random
import numpy as np
import argparse
from pysam import FastaFile
from processing import Error, FormatterIO, collect_args
from constants import Constants, Variant_Type
#from arguments import collect_args
#import tracemalloc   # only for testing
import sys
from enum import Enum
import time

# issues: insertion
# lowercase = invert
# _ = space
# ' = complement


time_start = time.time()

class Structural_Variant():
    def __init__(self, sv_type, lengths):
        '''
        sv_type: integer or tuple containing transformation (source, target)
        lengths: list containing tuple(s) (min_length, max_length)
        '''
        
        if isinstance(sv_type, int):
            self.type = Variant_Type(sv_type)
            self.source, self.target = Constants.SV_KEY[self.type]
            self.name = self.type.name
        else:
            self.type = None
            self.source, self.target = sv_type
            self.name = str(self.source) + "-" + str(self.target)

        # determine size of "letters" or pieces within variant
        self.lengths = []
        self.initialize_lengths(lengths)

        # initialize event classes
        self.events = []
        self.events_dict = dict()
        self.initialize_events()
        self.source_event_blocks = []
        self.target_event_blocks = []

        # in the event that sv is unable to be simulated due to random placement issues, will be turned on later
        self.active = False

        # 1 = homozygous, 0 = heterozygous
        self.ishomozygous = -1
        self.hap1 = -1
        self.hap2 = -1
    
    def __repr__(self):
        return "<SV transformation {} -> {} with lengths {}>".format(self.source, self.target, self.lengths)
    
    def initialize_lengths(self, lengths):
        if len(lengths) > 1:    # values given by user represents custom ranges for each event of variant
            assert (len(lengths) == len([letter for letter in self.source if letter != "-"])) 

            for leng in lengths:
                self.lengths.append(random.randint(leng[0], leng[1]))
            self.lengths = tuple(self.lengths)
            self.req_space = sum(self.lengths)

        elif len(lengths) == 1: # value given by user represents length (same range) of each event within variant
            count = 0
            num_char = len(self.source)
            for x in range(num_char):
                if self.source[x] != "-":
                    length = random.randint(lengths[0][0], lengths[0][1])
                else:
                    length = 0
                self.lengths.append(length)
                count += length
            
            self.req_space = count

        else:
            raise Exception("Lengths parameter expects at least one tuple")
        
    def initialize_events(self):
        for idx, letter in enumerate(self.source):
            event = Event(self, self.lengths[idx], letter)
            self.events.append(event)
            self.events_dict[letter] = event
        
        return self.events
    
    def generate_blocks(self):
        # groups together letters between dispersion events (_)
        # returns list of lists
        # Ex. AB_CD -> [["A","B"], ["C","D"]]
        def find_blocks(letter_transform):
            blocks = [[]]
            for letter in letter_transform:
                if letter != "_":
                    blocks[-1].append(letter)
                else:
                    blocks.append([])
            return blocks
        self.source_event_blocks = find_blocks(self.source)
        self.target_event_blocks = find_blocks(self.target)

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
        def generate_seq(length):
            # helper function for insertions
            rand_seq = ""
            for x in range(length):
                num = random.randint(1,4)
                if num == 1:
                    rand_seq += "A"
                elif num == 2:
                    rand_seq += "T"
                elif num == 3:
                    rand_seq += "G"
                elif num == 4:
                    rand_seq += "C"
            return rand_seq

        decode_dict = {"invert": lambda string: complement(string[::-1]),
                       "identity": lambda string: string,
                       "complement": complement} 
        encode = self.events_dict    # maps letter like A or B to base pairs on reference 
        #print("Encode_dict: ", encode)

        # find all blocks between dispersion events
        self.generate_blocks()
        
        changed_fragments = []
        for idx, block in enumerate(self.target_event_blocks):
            new_frag = ""
            for x, ele in enumerate(block):
                upper_str = ele.upper()         # changes all lowercase letters (if there are any) to uppercase so we can map to the nucleotide bases 
    
                if any(c.islower() for c in ele):   # checks if lowercase letters exist in ele
                    curr_piece = decode_dict["invert"](encode[upper_str[0]].source_frag)   # remember decode_dict stores functions as the values
    
                elif upper_str[0] in encode:
                    curr_piece = decode_dict["identity"](encode[upper_str[0]].source_frag)
                
                elif ele == "_":               # _ refers to a space
                    raise Exception("Dispersion event detected within block: {}".format(self.event_blocks))
                else:  # insertion case
                    curr_piece = generate_seq(20)
                
                new_frag += curr_piece
            start_event = encode[self.source_event_blocks[idx][0].upper()]
            end_event = encode[self.source_event_blocks[idx][-1].upper()]
            start_pos, end_pos = start_event.start, end_event.end
            changed_fragments.append([start_pos, end_pos, new_frag])
                
        #print("ref {} -> transformed {} for transformation {}".format(ref_piece, change_genome, self.type.value))
        self.changed_fragments = changed_fragments

        return changed_fragments

class Event():
    '''represents the letters/events within a SV transformation'''
    def __init__(self, sv_parent, length, letter):
        '''
        sv_parent: Structural Variant, event is a part of larger SV
        '''
        self.sv_parent = sv_parent
        self.length = length
        self.letter = letter # refers to letter in SV's transformation
        self.source_frag = ""
        self.target_frag = ""
        self.start = None
        self.end = None
    
    def __repr__(self):
        if self.start and self.end and self.source_frag:
            return "Event at {} with fragment {}".format((self.start, self.end), self.source_frag)
        elif self.start and self.end:
            return "Event at {} with parent {}".format((self.start, self.end), self.sv_parent)
        else:
            return "Event with parent {}".format(self.sv_parent)

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
        message = "SVs1: " + str([self.svs[x].type.name for x in range(len(self.svs1)) if self.svs1[x]]) + "\n"
        message += "SVs2: " + str([self.svs[y].type.name for y in range(len(self.svs2)) if self.svs2[y]]) + "\n"
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

    def produce_variant_genome(self, fasta1_out, fasta2_out, bedfile, initial_reset = True, verbose = False):
        '''
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        '''
        global time_start
        if initial_reset:
            self.formatter.reset_file(fasta1_out)
            self.formatter.reset_file(fasta2_out)

        ref_fasta = self.ref_fasta
        for id in self.order_ids:
            # edit chromosome
            self.rand_edit_svs(id, ref_fasta)
            print("Finished edits in {} seconds\n".format(time.time() - time_start))
            time_start = time.time()

            for x in range(2):
                if x == 0:
                    fasta_out = fasta1_out
                elif x == 1:
                    fasta_out = fasta2_out

                # account for homozygous and heterogeneous variants
                edits_x = []
                for sv in self.svs:
                    if sv.active and sv.hap[x]:
                        edits_x.extend(sv.changed_fragments)
                print("Edits_x: ", edits_x)
                
                # export edited chromosomes to FASTA files
                self.formatter.export_piece({id: edits_x}, fasta_out, ref_fasta, verbose = verbose)

                print("ID {} altered and saved in fasta file {} in {} seconds\n".format(id, fasta_out, time.time() - time_start))
                time_start = time.time()

            # export variant data to BED file
            active_svs = [sv for sv in self.svs if sv.active]
            ishomozygous = [sv.ishomozygous for sv in self.svs if sv.active]
            self.formatter.export_to_bedpe(active_svs, ishomozygous, id, bedfile, reset_file = initial_reset)

            initial_reset = False

        return True
    
    def rand_select_svs(self, svs, id, ref_fasta):
        def percent_N(seq):
            total = 0
            for char in seq:
                if char == "N":
                    total += 1
            return total / len(seq)
        def is_overlapping(event_ranges, addition, letter):
            # addition: tuple (start, end)
            # event_ranges: list containing tuples
            # checks if addition overlaps with any of the events already stored
            # NOTE: if addition is a dispersion event (letter == "_"), then automatically return False

            if letter != "_":
                for event in event_ranges:
                    if event[1] > addition[0] and event[0] < addition[1]:
                        return True
            return False
        
        self.event_ranges = []
        chr_len = self.len_dict[id]
        for sv in svs:
            #print("Current SV: ", sv)
            tries = 0 # number of attempts to place sv randomly
            valid = False
            while not valid:
                tries += 1
                valid = True
                if tries > 100 or chr_len - sv.req_space - 1 <= 0:
                    print("Failure to simulate SV \"{}\"".format(sv))
                    valid = False
                    break
                
                elif chr_len - sv.req_space - 1 > 0:
                    start_pos = random.randint(0, chr_len - sv.req_space - 1)
                    for sv_event in sv.events:

                        # check to see if chosen spot is a valid position
                        frag = ref_fasta.fetch(id, start_pos, start_pos + sv_event.length)
                        #print("Percent_N, {} -> {}".format(frag, percent_N(frag)))
                        #print("Is_overlapping, {} -> {}".format((self.event_ranges, (start_pos, start_pos + sv_event.length)), is_overlapping(self.event_ranges, (start_pos, start_pos + sv_event.length), sv_event.letter)))
                        if percent_N(frag) > 0.05 or is_overlapping(self.event_ranges, (start_pos, start_pos + sv_event.length), sv_event.letter):
                            valid = False

                        # store start and end position and reference fragment
                        sv_event.start, sv_event.end = start_pos, start_pos + sv_event.length
                        sv_event.source_frag = frag
                        start_pos += sv_event.length
            
            # dispersion events with a space are not off-limits to more events
            if valid:
                sv.active = True
                self.event_ranges.extend([(event.start, event.end) for event in sv.events if event.letter != "_"])
                sv.start = sv.events[0].start
                sv.end = sv.events[-1].end
            #print("\n")
                
        return self.event_ranges

    def rand_edit_svs(self, id, ref_fasta):
        # select random positions for SVs
        self.rand_select_svs(self.svs, id, ref_fasta)

        for sv in self.svs:
            if sv.active:
                # make edits and store in sv object
                sv.change_fragment()
        




if __name__ == "__main__":

    '''
    args = collect_args()
    #tracemalloc.start()
    fasta_in = args[0]
    yaml_in = args[1]
    fasta1_out = args[2]
    fasta2_out = args[3]
    bed_out = args[4]
    '''

    fasta_in = "debugging/inputs/test.fna"
    yaml_in = "debugging/inputs/par_test.yaml"
    #test_svs = [[12,[(30,100)]], [16,[(300,8000)]],[8,[(500,1000)]]]
    fasta1_out = "debugging/inputs/test1_out.fna"
    fasta2_out = "debugging/inputs/test2_out.fna"
    bed_out = "debugging/inputs/out.bed"

    sim = SV_Simulator(fasta_in, yaml_in)
    print(sim.svs)
    #print(sim.rand_edit_svs("Chromosome19", sim.ref_fasta))
    sim.produce_variant_genome(fasta1_out, fasta2_out, bed_out, verbose = False)
    #print("\n" + str(sim))

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    #tracemalloc.stop()


