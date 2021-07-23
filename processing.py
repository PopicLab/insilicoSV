import yaml
import sys
from constants import Constants, Operations
from pysam import FastaFile
#import tracemalloc

class Error(Exception):
    pass

class FormatterIO():
    def __init__(self, filein):
        self.bedpe_counter = 1
        self.bedpe_counter_dif = 1
        self.initialize_fasta_export(filein)

    def initialize_fasta_export(self, filein):

        # initialization for exporting to fasta files
        try:
            self.fin_export1 = open(filein, "r")   # needed to read original genome when exporting
            self.fin_export2 = open(filein, "r")   # needed to read original genome when exporting
        except FileNotFoundError:
            print("Error: File \"{}\" Not Found".format(filein))
        self.fout_export = False

    def yaml_to_var_list(self, filename):
        par_list = yaml.full_load(open(filename))
        variant_list = []
        for variant in par_list:
            var_ranges = []
            if isinstance(variant[Constants.MIN_LENGTH_ATTR], list):
                var_ranges = [(variant[Constants.MIN_LENGTH_ATTR][x], variant[Constants.MAX_LENGTH_ATTR][x]) for x in range(len(variant[Constants.MIN_LENGTH_ATTR]))]
            else:
                var_ranges = [(variant[Constants.MIN_LENGTH_ATTR], variant[Constants.MAX_LENGTH_ATTR])]
            variant_list.append([variant[Constants.TYPE_ATTR], variant[Constants.NUM_ATTR], var_ranges])
        return variant_list

    def find_lcs(self, str1, str2):
        # Finds longest common subsequence between str1 and str2 such that if there is a _ in str1 and str2, it MUST belong in the lcs
        # If there is a space (_) in str1 and str2, then simply find the lcs of before and after the _ and put them together
        # The lcs MUST have all the _ in str1 and str2 as a space cannot be transformed
        # str1 = source, str2 = target

        def lcs(str1, str2):
            # -> returns the longest common subsequence
            m, n = len(str1), len(str2)
        
            L = [[(0,"") for i in range(n + 1)]
                for i in range(m + 1)]
    
            # bottom-up approach
            for i in range(m + 1):
                for j in range(n + 1):
                    if (i == 0 or j == 0):
                        L[i][j] = (0,"")
                    elif(str1[i - 1].upper() == str2[j - 1].upper()):
                        L[i][j] = (L[i - 1][j - 1][0] + 1, L[i - 1][j - 1][1] + str2[j-1])
                    else:
                        L[i][j] = max(L[i - 1][j],
                                    L[i][j - 1])
        
            return L[m][n]

        # collects locations of _ (space) in str1 and str2
        loc1 = [-1]
        loc1.extend([index for index in range(len(str1)) if str1[index] == "_"])
        loc1.append(len(str1))

        loc2 = [-1]
        loc2.extend([index for index in range(len(str2)) if str2[index] == "_"])
        loc2.append(len(str2))

        if len(loc1) != len(loc2):
            raise Error("Unequal number of _ across str1 and str2, str1: {}, str2: {}".format(str1, str2))
        
        common = ""
        start = 0
        for end in range(1, len(loc1)):
            common += lcs(str1[loc1[start] + 1: loc1[end]], str2[loc2[start] + 1: loc2[end]])[1]
            if end != len(loc1) - 1:
                common += "_"
            
            start = end
        return common
    
    def export_to_bedpe(self, svs, ishomozygous, id, bedfile, reset_file = True):
        def write_to_file(sv, source_s, source_e, target_s, target_e, transform, ishomozygous, letter, order = 0):
            assert (letter != "_")

            # do not write transformations of size 0
            if source_e > source_s:
                with open(bedfile, "a") as fout:
                    row = [str(id),
                            str(source_s),
                            str(source_e + 1),
                            str(id),
                            str(target_s),
                            str(target_e + 1),
                            transform,
                            str(source_e - source_s),
                            str(int(ishomozygous)) + "/1",
                            sv.type.name,
                            str(self.bedpe_counter),
                            str(order)]

                    fout.write("\t".join(row) + "\n")

        if reset_file:
            self.reset_file(bedfile)
        
        for x in range(len(svs)):
            sv = svs[x]

            source, target = sv.source, sv.target
            letter_pos = dict()
            place = sv.start
            for y in range(len(sv.lengths)):
                letter_pos[source[y]] = (place, place + sv.lengths[y])
                place += sv.lengths[y]

            # longest common subsequence represents the letters that stay in place
            # the goal here is to minimize the number of transformations to get from source to target
            source, target = "".join(source), "".join(target)
            same_place = self.find_lcs(source, target)
            #print(source, target, same_place)

            # remove letters in source which do not remain in same place
            index = 0
            for sr_index in range(len(source)):
                if index >= len(same_place) or source[sr_index].upper() != same_place[index].upper():
                    start, end = letter_pos[source[sr_index].upper()]
                    write_to_file(sv, start, end, start, start, Operations.DEL.name, ishomozygous[x], source[sr_index])
                else:
                    index += 1
            
            # add letters from target not in same_place
            index = 0
            insert_point = sv.start
            order = 0
            for tr_index in range(len(target)):
                if index >= len(same_place) or target[tr_index].upper() != same_place[index].upper():
                    order += 1
                    start, end = letter_pos[target[tr_index].upper()]
                    write_to_file(sv, start, end, insert_point, insert_point, Operations.DUP.name, ishomozygous[x], target[tr_index], order)
                    
                    # Inversion applies in either case
                    if target[tr_index].islower():
                        write_to_file(sv, start, end, insert_point, insert_point, Operations.INV.name, ishomozygous[x], target[tr_index], order)
            
                elif target[tr_index].upper() == same_place[index].upper():
                    insert_point = letter_pos[same_place[index].upper()][1]
                    order = 0
                    index += 1
                
                    # Inversion applies in either case
                    if target[tr_index].islower():
                        start, end = letter_pos[target[tr_index].upper()]
                        write_to_file(sv, start, end, start, end, Operations.INV.name, ishomozygous[x], target[tr_index], order)
            
            self.bedpe_counter += 1
    
    def export_piece(self, variants, fasta_out, fasta_file, verbose = False):
        with open(fasta_out,"a") as fout_export:
            for id in variants:
                if id not in fasta_file.references:
                    raise KeyError("ID {} not found in inputted fasta file".format(id))
            
            for id in variants:
                if verbose:
                    print("New ID: ", id)
                fout_export.write(">" + str(id) + "\n")
                chr_variants = variants[id]
                chr_variants.sort()
                chr_variants.append([fasta_file.get_reference_length(id), fasta_file.get_reference_length(id), ""]) # to ensure that all bases after last variant are included

                pos = 0
                for variant in chr_variants:
                    var_start, var_end = variant[0], variant[1]
                    for idx in range(pos, var_end):
                        if idx < var_start:
                            c = fasta_file.fetch(id, idx, idx+1)
                            fout_export.write(c)

                        elif idx == var_start:
                            fout_export.write(variant[2])

                        if idx == var_end - 1 and var_start == var_end:   # for insertions, insert piece right before position var_end
                            fout_export.write(variant[2])
                    pos = var_end
                fout_export.write("\n")

    def reset_file(self, filename):
        #print("Overwritting File {}...".format(filename))
        with open(filename, "w") as f_reset:
            f_reset.truncate()


if __name__ == "__main__":
    pass
    
    #tracemalloc.start()
    
    '''x = 10
    y = 10
    filein = "debugging/inputs/test.fna"
    fileout1 = "debugging/inputs/test1_tmp_out.fna"
    fileout2 = "debugging/inputs/test2_tmp_out.fna"
    fasta = FormatterIO(filein)
    fasta_file = FastaFile(filein)


    order_ids = fasta_file.references
    len_dict = dict()
    for id in order_ids:
        len_dict[id] = fasta_file.get_reference_length(id)
    print("Len_dict: ", len_dict)
    #fasta.next()
    #print(fasta.fetch(80,85))
    variants = {"Chromosome19": [[5,5,""]]}
    FastaFile.reset_file(fileout1)
    fasta.export_piece(variants, fileout1, 0, verbose = False)
    variants = {"Chromosome19": [[2,5,"UUUUU"],[7,10, "U"], [15, 20, "TC"], [75,77,"UUUUUU"]]}
    variants2 = {"Chromosome19": [[2,5,"UUUUUU"]], "Chromosome21": [[4,15,"AAA"]]}
    fasta.reset_file(fileout1)
    fasta.reset_file(fileout2)
    fasta.export_piece_fasta(variants, fileout1, fasta_file, verbose = True)
    fasta.export_piece_fasta(variants2, fileout2, fasta_file, verbose = True)
    #fasta.export_piece(variants, fileout1, 0, len_dict, verbose = True)

    print(fasta_file.fetch("Chromosome19", 1, 2))
    #fasta.close()'''
    # CTUUUUUTCUCTAGATCCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCT


    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    #tracemalloc.stop()

    '''tracemalloc.start()
    with open("test.txt", "r") as fin:
        a = fin.readline()
        fin.readline()

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()'''


