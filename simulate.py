import random
import numpy as np
import argparse
from pysam import FastaFile
from processing import Error, FormatterIO
from constants import Constants, Variant_Type
from arguments import collect_args
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
    def __init__(self, sv_type, lengths, ishomozygous = None):
        '''
        sv_type: integer
        lengths: list containing tuple(s) (min_length, max_length)
        '''
        
        self.type = Variant_Type(sv_type)
        self.source, self.target = Constants.SV_KEY[self.type]
        self.lengths = []
        # self.piece_des = bedpe_key[self.type]
        #self.ishomozygous = ishomozygous

        # determine size of "letters" or pieces within variant
        if len(lengths) > 1:    # values given by user represents length of each piece of variant
            assert (len(lengths) == len(Constants.SV_KEY[self.type][0])) 

            for leng in lengths:
                self.lengths.append(random.randint(leng[0], leng[1]))
            self.lengths = tuple(self.lengths)
            self.req_space = sum(self.lengths)

        elif len(lengths) == 1: # value given by user represents length of each event within variant
            count = 0
            num_char = len(Constants.SV_KEY[self.type][0])
            for x in range(num_char):
                length = random.randint(lengths[0][0], lengths[0][1])
                self.lengths.append(length)
                count += length
            
            self.req_space = count

        else:
            raise ValueError("Lengths parameter expects at least one tuple")
    
    def __str__(self):
        return "SV Object with Type {} and Variant Lengths {}".format(self.type, self.lengths)

    def decode(self, begin, end, ref_piece):
        '''
        begin: tuple
        end: tuple, shows ending transformation from variant
        '''
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
        encode = dict()    # maps letter like A or B to base pairs on reference 
        self.target_lengths = []
        
        assert (len(begin) == len(self.lengths))

        curr_index = 0
        for x in range(len(begin)):
            ele = begin[x]
            encode[ele] = ref_piece[curr_index: curr_index + self.lengths[x]]
            curr_index += self.lengths[x]
        
        #print("Encode_dict: ", encode)

        change_genome = ""
        for y in range(len(end)):
            ele = end[y]
            upper_str = ele.upper()         # changes all lowercase letters (if there are any) to uppercase so we can map to the nucleotide bases 

            if any(c.islower() for c in ele):   # checks if lowercase letters exist in ele
                curr_piece = decode_dict["invert"](encode[upper_str[0]])   # remember decode_dict stores functions as the values

            elif ele == "_":               # _ refers to a space
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            else:
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            self.target_lengths.append(len(curr_piece))
            change_genome += curr_piece
        #print("ref {} -> transformed {} for transformation {}".format(ref_piece, change_genome, self.type.value))

        return change_genome

    def change_piece(self, ref):
        return self.decode(self.source, self.target, ref)
        #print("SV of Type {} starting at index {} and ending at index {}".format(self.type, self.start, self.end))


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
        svs = self.formatter.yaml_to_var_list(par_file)
        self.initialize_svs(svs)

        print("Finished Setting up Simulator in {} seconds\n".format(time.time() - time_start))
        time_start = time.time()
    
    def __repr__(self):
        message = "SVs1: " + str([self.svs[x].type.name for x in range(len(self.svs1)) if self.svs1[x]]) + "\n"
        message += "SVs2: " + str([self.svs[y].type.name for y in range(len(self.svs2)) if self.svs2[y]]) + "\n"
        return message
    
    def initialize_svs(self, svs):
        self.svs = []
        self.svs1 = []
        self.svs2 = []
        for sv in svs:
            for num in range(sv[1]):
                self.svs.append(Structural_Variant(sv[0], sv[2])) # inputs: SV type, range of lengths
                draw = random.randint(1,4)
                if draw == 3 or draw == 4:   # sv applies to both haplotypes
                    self.svs1.append(1)
                    self.svs2.append(1)
                elif draw == 2:
                    self.svs1.append(1)
                    self.svs2.append(0)
                elif draw == 1:
                    self.svs1.append(0)
                    self.svs2.append(1)
        self.ishomozygous = [self.svs1[index] == 1 and self.svs2[index] == 1 for index in range(len(self.svs1))]

    def export_variant_genome(self, fasta1_out, fasta2_out, bedfile, initial_reset = True, verbose = False):
        '''
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        '''
        global time_start
        if initial_reset:
            self.formatter.reset_file(fasta1_out)
            self.formatter.reset_file(fasta2_out)

        ref_fasta = self.ref_fasta
        for id in self.order_ids:
            # shuffle svs, now assume the SVs will be ordered in this way in the altered chromosome
            temp = list(zip(self.svs, self.svs1, self.svs2))
            random.shuffle(temp)
            self.svs, self.svs1, self.svs2 = zip(*temp)

            # edit chromosome
            edits = self.rand_edit_svs(self.svs, id, ref_fasta)
            print("Finished edits in {} seconds\n".format(time.time() - time_start))
            time_start = time.time()

            for x in range(2):
                if x == 0:
                    fasta_out = fasta1_out
                    activations = self.svs1
                elif x == 1:
                    fasta_out = fasta2_out
                    activations = self.svs2

                # account for homozygous and heterogeneous variants
                edits_x = [edits[val] for val in range(len(activations)) if activations[val]]
                
                # export edited chromosomes to FASTA files
                self.formatter.export_piece({id: edits_x}, fasta_out, ref_fasta, verbose = verbose)

                print("ID {} altered and saved in fasta file {} in {} seconds\n".format(id, fasta_out, time.time() - time_start))
                time_start = time.time()

            # export variant data to BED file
            self.formatter.export_to_bedpe(self.svs, self.ishomozygous, id, bedfile, reset_file = initial_reset)

            initial_reset = False

        return True


    def rand_edit_svs(self, svs, id, ref_fasta):
        '''
        -> tuples with random start and end positions referring to corresponding sv from input list (including start, excluding end)
        '''
        sum_spaces = 0
        for sv in svs:
            sum_spaces += sv.req_space

        avail_spaces = self.len_dict[id] - sum_spaces - 1
        if avail_spaces < 0:
            raise ValueError("Size of variants is too big for the chromosome id {}!".format(id))

        curr_place = 0
        edited_pieces = []

        for x in range(len(svs)):             # think of this as randomly inserting spaces in between the ordered SVs
            # determine spacing between variants
            curr_sv = svs[x]
            center, sd = avail_spaces / (len(svs) - x + 1), avail_spaces / 8
            space = int(np.clip(np.round(np.random.normal(center, sd,1)), 2, avail_spaces - 2))
            curr_place += space

            # edit the appropriate piece of the genome
            ref_piece = ref_fasta.fetch(id, curr_place, curr_place + curr_sv.req_space)
            edit = curr_sv.change_piece(ref_piece)

            # saving info to export into fasta and bed file
            curr_sv.start = curr_place
            curr_sv.end = curr_place + curr_sv.req_space
            curr_sv.len_edit = len(edit)
            #edit = "U" * curr_sv.len_edit    # ONLY FOR TESTING AND DEVELOPMENT
            edited_pieces.append([curr_place, curr_place + curr_sv.req_space, edit])

            curr_place += curr_sv.req_space
            avail_spaces -= space

        return edited_pieces
    

if __name__ == "__main__":


    args = collect_args()
    #tracemalloc.start()
    fasta_in = args[0]
    yaml_in = args[1]
    fasta1_out = args[2]
    fasta2_out = args[3]
    bed_out = args[4]

    '''fasta_in = "debugging/inputs/test.fna"
    yaml_in = "par_test.yaml"
    #test_svs = [[12,[(30,100)]], [16,[(300,8000)]],[8,[(500,1000)]]]
    fasta1_out = "debugging/inputs/test1_out.fna"
    fasta2_out = "debugging/inputs/test2_out.fna"
    bed_out = "debugging/inputs/out.bed"'''

    sim = SV_Simulator(fasta_in, yaml_in)
    sim.export_variant_genome(fasta1_out, fasta2_out, bed_out, verbose = False)
    #print("\n" + str(sim))

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    #tracemalloc.stop()



