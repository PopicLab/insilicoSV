import random
import numpy as np
from processing import Error, FastaFile, Formater
import tracemalloc   # only for testing
import sys
from enum import Enum

# only set for testing and development
#random.seed(10)
#np.random.seed(10)

class Variant_Type(Enum):
    INS = 1
    DEL = 2
    INV = 3
    DUP = 4
    TRANS = 5
    dupINVdup = 6
    delINVdel = 7
    delINVdup = 8
    dupINVdel = 9
    delINV = 10
    INVdel = 11
    dDUP_iDEL = 12
    INS_iDEL = 13
    dupINV = 14
    INVdup = 15
    dDUP = 16

bedpe_key = {Variant_Type.INS: ["INS"],     
            Variant_Type.DEL: ["DEL"],
            Variant_Type.INV: ["INV"],
            Variant_Type.DUP: ["DUP"],
            Variant_Type.TRANS: ["TRANS", "TRANS", "TRANS"],
            Variant_Type.dupINVdup: ["DUP", "INV", "DUP"],
            Variant_Type.delINVdel: ["DEL", "INV", "DEL"],
            Variant_Type.delINVdup: ["DEL", "INV", "DUP"],
            Variant_Type.dupINVdel: ["DUP", "INV", "DEL"],
            Variant_Type.delINV: ["DEL", "INV"],
            Variant_Type.INVdel: ["INV", "DEL"],
            Variant_Type.dDUP_iDEL: ["DUP", "d", "DEL"],
            Variant_Type.INS_iDEL: ["INS", "DEL"],
            Variant_Type.dupINV: ["DUP", "INV"],
            Variant_Type.INVdup: ["INV", "DUP"],
            Variant_Type.dDUP: ["DUP", "d"]}

# issues: insertion
# lowercase = invert
# _ = space
# ' = complement

sv_key = {Variant_Type.INS: [("A",), ("A",)],     
            Variant_Type.DEL: [("A",), ()],
            Variant_Type.INV: [("A",), ("a",)],
            Variant_Type.DUP: [("A",), ("A","A")],
            Variant_Type.TRANS: [("A","_","B"), ("B","_","A")],
            Variant_Type.dupINVdup: [("A","B","C"), ("A","c","b","a","C")],
            Variant_Type.delINVdel: [("A","B","C"), ("b",)],
            Variant_Type.delINVdup: [("A","B","C"), ("c","b","C")],
            Variant_Type.dupINVdel: [("A","B","C"), ("A","b","a")],
            Variant_Type.delINV: [("A","B"), ("b",)],
            Variant_Type.INVdel: [("A","B"), ("a",)],
            Variant_Type.dDUP_iDEL: [("A","_","B"), ("A","_","A")],
            Variant_Type.INS_iDEL: [("A","_","B"), ("_","A")],
            Variant_Type.dupINV: [("A","B"), ("A","b","a")],
            Variant_Type.INVdup: [("A","B"), ("b'","a","B")],
            Variant_Type.dDUP: [("A","_"), ("A","_","A")]}


class Structural_Variant():
    def __init__(self, sv_type, lengths, ishomogeneous = None):
        '''
        sv_type: integer
        lengths: list containing tuple(s) (min_length, max_length)
        '''
        self.type = Variant_Type(sv_type)
        self.lengths = []
        self.piece_des = bedpe_key[self.type]
        #self.ishomogeneous = ishomogeneous

        # determine size of "letters" or pieces within variant
        if len(lengths) > 1:    # values given by user represents length of each piece of variant
            assert (len(lengths) == len(sv_key[self.type][0])) 

            for leng in lengths:
                self.lengths.append(random.randint(leng[0], leng[1]))
            self.lengths = tuple(self.lengths)
            self.req_space = sum(self.lengths)

        elif len(lengths) == 1: # value given by user represents total length of variant
            total_length = random.randint(lengths[0][0], lengths[0][1])
            self.req_space = total_length
            count = 0
            num_char = len(sv_key[self.type][0])
            for x in range(num_char):
                if x == num_char - 1:
                    self.lengths.append(total_length - count)
                else:
                    self.lengths.append(total_length // num_char)
                    count += total_length // num_char

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
            #print("Ele: ", ele)
            upper_str = ele.upper()         # changes all lowercase letters (if there are any) to uppercase so we can map to the nucleotide bases 

            if any(c.islower() for c in ele):   # checks if lowercase letters exist in ele
                curr_piece = decode_dict["invert"](encode[upper_str[0]])   # remember decode_dict stores functions as the values

            elif ele == "_":               # _ refers to a space
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            else:
                curr_piece = decode_dict["identity"](encode[upper_str[0]])
            
            #if "'" in ele:
            #    curr_piece = decode_dict["complement"](curr_piece)
            self.target_lengths.append(len(curr_piece))
            change_genome += curr_piece
        #print("ref {} -> transformed {} for transformation {}".format(ref_piece, change_genome, self.type.value))

        return change_genome

    def change_piece(self, ref):
        return self.decode(sv_key[self.type][0], sv_key[self.type][1], ref)
        #print("SV of Type {} starting at index {} and ending at index {}".format(self.type, self.start, self.end))


class SV_Simulator():
    def __init__(self, ref_file, par_file):
        '''
        ref: fasta filename to reference genome
        svs: list (type, number, range of length of variant)

        '''

        print("Setting Up Simulator...")
        self.ref_file = ref_file
        self.ref_fasta = FastaFile(ref_file) 
        print("Length Dict: ", self.ref_fasta.len_dict)
        self.formater = Formater()
        svs = self.formater.yaml_to_var_list(par_file)
        random.shuffle(svs)     # now assume the SVs will be ordered in this way in the altered chromosome

        self.svs = []
        self.svs1 = []
        self.svs2 = []
        print("Setting Up Structural Variants")
        for sv in svs:
            for num in range(sv[1]):
                self.svs.append(Structural_Variant(sv[0], sv[2]))
                draw = random.randint(1,3)
                if draw == 3:   # sv applies to both haplotypes
                    self.svs1.append(1)
                    self.svs2.append(1)
                elif draw == 2:
                    self.svs1.append(1)
                    self.svs2.append(0)
                elif draw == 1:
                    self.svs1.append(0)
                    self.svs2.append(1)

        print("Finished Setting up Simulator")

    
    def __str__(self):
        message = "SVs1: " + str([self.svs[x].type.name for x in range(len(self.svs1)) if self.svs1[x]]) + "\n"
        message += "SVs2: " + str([self.svs[y].type.name for y in range(len(self.svs2)) if self.svs2[y]]) + "\n"
        return message

    def export_variant_genome(self, fasta1_out, fasta2_out, bedfile_list, initial_reset = True, verbose = False):
        '''
        bedfile_list: list (if it contains more than one element, export to the different bed files)
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        '''

        if initial_reset:
            FastaFile.reset_file(fasta1_out)
            FastaFile.reset_file(fasta2_out)

        ref_fasta = self.ref_fasta
        for id in self.ref_fasta.order_ids:
            # reshuffle svs
            temp = list(zip(self.svs, self.svs1, self.svs2))
            random.shuffle(temp)
            self.svs, self.svs1, self.svs2 = zip(*temp)

            # edit chromosome
            edits = self.rand_select_svs(self.svs, id, ref_fasta)

            # move on to next chromosome id
            ref_fasta.next()

            for x in range(2):
                bedfile = bedfile_list[0]
                if x == 0:
                    fasta_out = fasta1_out
                    activations = self.svs1
                elif x == 1:
                    fasta_out = fasta2_out
                    activations = self.svs2

                if len(bedfile_list) > 1:
                    bedfile = bedfile_list[x]

                # account for homogeneous and heterogeneous variants
                edits_x = [edits[val] for val in range(len(activations)) if activations[val]]
                
                # export edited chromosomes to FASTA files
                ref_fasta.export_piece({id: edits_x}, fasta_out, x, verbose = verbose)

                print("ID {} altered and saved in fasta file {}".format(id, fasta_out))

            # export variant data to BED file
            #print("SVS1: {}".format(self.svs1))
            #print("SVS2: {}".format(self.svs2))
            ishomogeneous = [self.svs1[index] == 1 and self.svs2[index] == 1 for index in range(len(self.svs1))]
            self.formater.export_to_bedpe(self.svs, ishomogeneous, id, bedfile, reset_file = initial_reset)
            #self.formater.export_to_bed12(svs, id, self.ref_file, fasta_out, bedfile, reset_file = initial_reset)

            initial_reset = False

        return True


    def rand_select_svs(self, svs, id, ref_fasta):
        '''
        -> tuples with random start and end positions referring to corresponding sv from input list (including start, excluding end)
        '''
        sum_spaces = 0
        for sv in svs:
            sum_spaces += sv.req_space

        avail_spaces = ref_fasta.len_dict[id] - sum_spaces - 1
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
            ref_piece = ref_fasta.fetch(curr_place, curr_place + curr_sv.req_space)
            edit = curr_sv.change_piece(ref_piece)

            # saving info to export into fasta and bed file
            curr_sv.start = curr_place
            curr_sv.end = curr_place + curr_sv.req_space
            curr_sv.len_edit = len(edit)
            edit = "U" * curr_sv.len_edit    # ONLY FOR TESTING AND DEVELOPMENT
            edited_pieces.append([curr_place, curr_place + curr_sv.req_space, edit])

            curr_place += curr_sv.req_space
            avail_spaces -= space
        #print("Size of SV: {} MB".format(sys.getsizeof(curr_sv)/10**6))

        return edited_pieces
    
    def determine_target_pos(self, svs):
        '''
        Determines position of edited pieces on altered chromosome given list of Stuctural Variants
        '''
        svs.sort(key = lambda sv: sv.start)
        offset = 0

        for sv in svs:
            # original chromosome blocks
            block_starts = []
            place = sv.start
            for x in range(len(sv.lengths)):
                block_starts.append(str(int(place)))
                place += sv.lengths[x]
            sv.block_starts = ','.join(block_starts)

            # edited chromosome
            sv.target_start = sv.start + offset
            sv.target_end = sv.target_start + sv.len_edit
            place = sv.target_start
            target_block_starts = []
            for y in range(len(sv.target_lengths)):
                target_block_starts.append(str(int(place)))
                place += sv.target_lengths[y]
            sv.target_block_starts = ','.join(target_block_starts)

            offset = sv.target_end - sv.end





if __name__ == "__main__":
    tracemalloc.start()

    fasta_in = "reference/test.fna"
    yaml_in = "par.yaml"
    #test_svs = [[12,[(30,100)]], [16,[(300,8000)]],[8,[(500,1000)]]]
    fasta1_out = "reference/test1_out.fna"
    fasta2_out = "reference/test2_out.fna"
    bed_out = "reference/out.bed"
    sim = SV_Simulator(fasta_in, yaml_in)
    sim.export_variant_genome(fasta1_out, fasta2_out, [bed_out], verbose = False)
    print("\n" + sim)
    #print(len(sim.svs1))
    #print(sim.svs1[5].__dict__)

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()



