import yaml
import sys
from constants import Constants, Operations
from pysam import FastaFile
import argparse
import os
#import tracemalloc

class Error(Exception):
    pass

class Config():
    def __init__(self):
        self.variant_config_list = []
    
    def run_checks(self, config):
        if Constants.TYPE_ATTR in config and Constants.TRANSFORM_SOURCE_ATTR in config and Constants.TRANSFORM_TARGET_ATTR in config:
            raise Exception("Only a type attribute (integer) or a source and target attribute can be present for SV, not both")
        elif Constants.TYPE_ATTR in config and not isinstance(config[Constants.TYPE_ATTR], int):
            raise Exception("Invalid {} type for SV \'type\' attribute, int expected".format(type(config[Constants.TYPE_ATTR])))
    
    def create_entries(self, config_entries):
        variant_config_list = []
        for variant_config in config_entries:
            self.run_checks(variant_config)

            # handles length of variant_config
            var_ranges = []
            if isinstance(variant_config[Constants.MIN_LENGTH_ATTR], list):
                var_ranges = [(variant_config[Constants.MIN_LENGTH_ATTR][x], variant_config[Constants.MAX_LENGTH_ATTR][x]) for x in range(len(variant_config[Constants.MIN_LENGTH_ATTR]))]
            else:
                var_ranges = [(variant_config[Constants.MIN_LENGTH_ATTR], variant_config[Constants.MAX_LENGTH_ATTR])]
            
            # type of variant_config (including custom)
            if Constants.TYPE_ATTR in variant_config:
                variant_config_list.append([variant_config[Constants.TYPE_ATTR], variant_config[Constants.NUM_ATTR], var_ranges])
            elif Constants.TRANSFORM_SOURCE_ATTR in variant_config and Constants.TRANSFORM_TARGET_ATTR in variant_config:
                variant_config_list.append([[tuple(variant_config[Constants.TRANSFORM_SOURCE_ATTR]), tuple(variant_config[Constants.TRANSFORM_TARGET_ATTR])], variant_config[Constants.NUM_ATTR], var_ranges])
        
        self.variant_config_list = variant_config_list
        return variant_config_list

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
        # config_entries = yaml.load(file, Loader=yaml.FullLoader)
        try:
            config_entries = yaml.full_load(open(filename))
        except:
            raise Exception("YAML File {} failed to be open".format(filename))

        config = Config()
        sv_configs = config.create_entries(config_entries)

        return sv_configs

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
            raise Exception("Unequal number of _ across str1 and str2, str1: {}, str2: {}".format(str1, str2))
        
        common = ""
        start = 0
        for end in range(1, len(loc1)):
            common += lcs(str1[loc1[start] + 1: loc1[end]], str2[loc2[start] + 1: loc2[end]])[1]
            if end != len(loc1) - 1:
                common += "_"
            
            start = end
        return common
    
    def export_to_bedpe(self, svs, id, bedfile, reset_file = True):
        def write_to_file(sv, source_s, source_e, target_s, target_e, transform, symbol, order = 0):
            assert (not symbol.startswith("_"))

            # do not write transformations of size 0
            if source_e > source_s or source_e == -1:
                with open(bedfile, "a") as fout:
                    row = [str(id),
                            str(source_s),
                            str(source_e + 1),
                            str(id),
                            str(target_s),
                            str(target_e + 1),
                            transform,
                            str(source_e - source_s),
                            str(int(sv.ishomozygous)) + "/1",
                            sv.name,
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
            for y, event in enumerate(sv.source_events):
                letter_pos[source[y]] = (place, place + event.length)
                place += event.length

            # longest common subsequence represents the letters that stay in place
            # the goal here is to minimize the number of transformations to get from source to target
            source, target = "".join(source), "".join(target)
            same_place = self.find_lcs(source, target)
            #print(source, target, same_place)

            # remove letters in source which do not remain in same place
            index = 0
            for sr_symbol in source:
                if index >= len(same_place) or sr_symbol.upper() != same_place[index].upper():
                    start, end = letter_pos[sr_symbol.upper()]
                    write_to_file(sv, start, end, start, start, Operations.DEL.name, sr_symbol)
                else:
                    index += 1
            
            # add letters from target not in same_place
            index = 0
            insert_point = sv.start
            order = 0
            for tr_symbol in target:
                if index >= len(same_place) or tr_symbol.upper() != same_place[index].upper():
                    order += 1
                    if tr_symbol.upper() in letter_pos:  # duplication
                        start, end = letter_pos[tr_symbol.upper()]
                        event_name = Operations.DUP.name
                    else:  # insertion detected
                        start, end = -1, -1
                        event_name = Operations.INS.name
                    write_to_file(sv, start, end, insert_point, insert_point, event_name, tr_symbol, order)
                    
                    # Inversion applies in either case
                    if tr_symbol.islower():
                        write_to_file(sv, start, end, insert_point, insert_point, Operations.INV.name, tr_symbol, order)
            
                elif tr_symbol.upper() == same_place[index].upper():
                    insert_point = letter_pos[same_place[index].upper()][1]
                    order = 0
                    index += 1
                
                    # Inversion applies in either case
                    if tr_symbol.islower():
                        start, end = letter_pos[tr_symbol.upper()]
                        write_to_file(sv, start, end, start, end, Operations.INV.name, tr_symbol, order)
            
            self.bedpe_counter += 1
    
    def export_piece(self, id, edits, fasta_out, fasta_file, verbose = False):
        with open(fasta_out,"a") as fout_export:
            if id not in fasta_file.references:
                raise KeyError("ID {} not found in inputted fasta file".format(id))
            
            if verbose:
                print("New ID: ", id)
            fout_export.write(">" + str(id) + "\n")
            chr_variants = list(edits)
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

def collect_args():
    parser = argparse.ArgumentParser(description='insilicoSV is a software to design and simulate complex structural variants, both novel and known.')
    parser.add_argument("ref", help="FASTA reference file")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("hap1", help="First output FASTA file (first haplotype)")
    parser.add_argument("hap2", help = "Second output FASTA file (second haplotype)")
    parser.add_argument("bedpe", help = "BEDPE file location to store variant info")
    parser.add_argument("-r", "--root", action="store", metavar="DIR", dest="root_dir", help="root directory for all files given as positional arguments")

    args = parser.parse_args()
    io = [args.ref, args.config, args.hap1, args.hap2, args.bedpe]
    if args.root_dir:
        io = [os.path.join(args.root_dir, ele) for ele in io]
    return io

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


