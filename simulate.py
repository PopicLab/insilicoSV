import random
import numpy as np
#from colorama import init
#from colorama import Fore
from processing import FastaFile, Formater
import tracemalloc   # only for testing
import sys

random.seed(10)

'''
1 = Insertion
2 = Deletion
3 = Inversion
4 = Duplication
5 = Translocation
6 = dupINVdup
7 = delINVdel
8 = delINVdup
9 = dupINVdel
10 = delINV
11 = INVdel
12 = dDUP-iDEL
13 = INS-iDEL
14 = dupINV
15 = INVdup
16 = dDUP
'''


class Structural_Variant():
    def __init__(self, sv_type, lengths):
        '''
        sv_type: integer
        lengths: list containing tuple(s) (min_length, max_length)
        '''

        # issues: insertion
        # lowercase = invert
        # _ = space
        # ' = complement

        self.sv_key_int = {1: [("A",), ("A",)],     
                           2: [("A",), ("",)],
                           3: [("A",), ("a'",)],
                           4: [("A",), ("A","A")],
                           5: [("A","_","B"), ("B","_","A")],
                           6: [("A","B","C"), ("A","c'","b","a'","C")],
                           7: [("A","B","C"), ("b",)],
                           8: [("A","B","C"), ("c'","b","C")],
                           9: [("A","B","C"), ("A","b","a'")],
                           10: [("A","B"), ("b",)],
                           11: [("A","B"), ("a",)],
                           12: [("A","_","B"), ("A","_","a'")],
                           13: [("A","_","B"), ("_","A")],
                           14: [("A","B"), ("A","b","a'")],
                           15: [("A","B"), ("b'","a","B")],
                           16: [("A","_"), ("A","_","a'")]}

        self.type = sv_type
        self.length = []
        if len(lengths) > 1:    # values given by user represents length of each piece of variant
            assert (len(lengths) == len(self.sv_key_int[sv_type][0]))

            for leng in lengths:
                self.length.append(random.randint(leng[0], leng[1]))
            self.length = tuple(self.length)
            self.req_space = sum(self.length)

        elif len(lengths) == 1: # value given by user represents total length of variant
            total_length = random.randint(lengths[0][0], lengths[0][1])
            self.req_space = total_length
            count = 0
            num_char = len(self.sv_key_int[sv_type][0])
            for x in range(num_char):
                if x == num_char - 1:
                    self.length.append(total_length - count)
                else:
                    self.length.append(total_length // num_char)
                    count += total_length // num_char

        else:
            raise ValueError("Lengths parameter expects at least one tuple")

    
    def __str__(self):
        return "SV Object with Type {} and Variant Lengths {}".format(self.type, self.length)

    def decode(self, begin, end, ref_piece):
        '''
        begin: tuple
        end: tuple, shows ending transformation from variant
        '''
        def complement(bases):
            output = ""
            for base in bases:
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
            
            return output

        decode_dict = {"invert": lambda string: string[::-1],
                       "identity": lambda string: string,
                       "complement": complement} 
        encode = dict()    # maps letter like A or B to base pairs on reference 
        
        assert (len(begin) == len(self.length))

        curr_index = 0
        for x in range(len(begin)):
            ele = begin[x]
            encode[ele] = ref_piece[curr_index: curr_index + self.length[x]]
            curr_index += self.length[x]
        
        #print("Encode_dict: ", encode)

        change_genome = ""
        for y in range(len(end)):
            ele = end[y]
            #print("Ele: ", ele)
            upper_str = ele.upper()         # changes all lowercase letters (if there are any) to uppercase so we can map to the nucleotide bases 

            if any(c.islower() for c in ele):   # checks if lowercase letters exist in ele
                curr_piece = decode_dict["invert"](encode[upper_str[0]])   # remember decode_dict stores functions as the values

            elif ele == "_":               # _ refers to a space
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            else:
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            if "'" in ele:
                curr_piece = decode_dict["complement"](curr_piece)
            change_genome += curr_piece
        #print("ref {} -> transformed {}".format(ref_piece, change_genome))

        return change_genome

    def change_piece(self, ref):
        return self.decode(self.sv_key_int[self.type][0], self.sv_key_int[self.type][1], ref)
        #print("SV of Type {} starting at index {} and ending at index {}".format(self.type, self.start, self.end))


class SV_Simulator():
    def __init__(self, ref_file, par_file):
        '''
        ref: fasta filename to reference genome
        svs: list (type, length of variant)

        *If there are multiple components to the variant, then the length parameter will be a tuple instead of an int
        '''
        self.ref_file = ref_file
        self.ref_fasta = FastaFile(ref_file)
        self.formater = Formater()
        svs = self.formater.yaml_to_dict(par_file)

        self.svs = []
        random.shuffle(svs)     # now assume the SVs will be ordered in this way in the altered genome
        for sv in svs:
            self.svs.append(Structural_Variant(sv[0], sv[1]))
    
    def export_variant_genome(self, fasta_out, bedfile = None):

        print("Length Dict: ", self.ref_fasta.len_dict)
        initial_reset = True
        for id in self.ref_fasta.len_dict:

            edits = self.rand_select_svs(self.ref_fasta.len_dict[id])

            if initial_reset:
                self.ref_fasta.export_piece({id: edits}, fasta_out)
                initial_reset = False
            else:
                self.ref_fasta.export_piece({id: edits}, fasta_out, reset_file=False)

            self.ref_fasta.next()   # move on to the next chromosome
            print("ID {} altered and saved".format(id))

        return True


    def rand_select_svs(self, ref_length):
        '''
        -> tuples with random start and end positions referring to corresponding sv from input list (including start, excluding end)
        '''
        sum_spaces = 0
        for sv in self.svs:
            sum_spaces += sv.req_space

        avail_spaces = ref_length - sum_spaces - 1
        if avail_spaces < 0:
            raise ValueError("Size of variants is too big for the given genome!")

        curr_place = 0
        edited_pieces = []
        for x in range(len(self.svs)):             # think of this as randomly inserting spaces in between the ordered SVs
            # determine spacing between variants
            curr_sv = self.svs[x]
            center, sd = avail_spaces / (len(self.svs) - x + 1), avail_spaces / 8
            space = float(np.clip(np.round(np.random.normal(center, sd,1)), 2, avail_spaces - 2))
            curr_place += space

            # edit the appropriate piece of the genome
            ref_piece = self.ref_fasta.fetch(curr_place, curr_place + curr_sv.req_space)
            edit = curr_sv.change_piece(ref_piece)
            edited_pieces.append([curr_place, curr_place + curr_sv.req_space, edit])

            curr_place += curr_sv.req_space
            avail_spaces -= space
        return edited_pieces


if __name__ == "__main__":
    tracemalloc.start()

    fasta_in = "test.txt"
    yaml_in = "par.yaml"
    #test_svs = [[12,[(30,100)]], [16,[(300,8000)]],[8,[(500,1000)]]]
    fasta_out = "test_out.txt"
    sim = SV_Simulator(fasta_in, yaml_in)
    sim.export_variant_genome(fasta_out)

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()



