import yaml
import sys
from constants import *
from pysam import FastaFile
import argparse
import os

class Config():
    def __init__(self,**entries):
        self.__dict__.update(entries)
        self.run_checks()
        if "fail_if_placement_issues" not in self.__dict__:
            self.fail_if_placement_issues = False
            #self.__dict__.update({"fail_if_placement_issues": False})
            #print(self.__dict__)

        for sv in self.SVs:
            # handles cases where user enters length range for all components within SV or specifies different ranges
            if isinstance(sv[MIN_LENGTH_ATTR], int):
                sv[RANGES_ATTR] = [(sv[MIN_LENGTH_ATTR], sv[MAX_LENGTH_ATTR])]
            else:
                sv[RANGES_ATTR] = list(zip(sv[MIN_LENGTH_ATTR], sv[MAX_LENGTH_ATTR]))
            
            # interchromosomal transformations (can only take place across dispersions)
            if "allow_interchromosomal" not in sv:
                sv["allow_interchromosomal"] = False
            
            # supply filler values in for source and target if the type is custom
            sv[TYPE_ATTR] = Variant_Type(sv[TYPE_ATTR])
            if sv[TYPE_ATTR] != Variant_Type.Custom:
                sv[TRANSFORM_SOURCE_ATTR] = None
                sv[TRANSFORM_TARGET_ATTR] = None

    def run_checks(self):
        config_svs = self.SVs
        for config_sv in config_svs:
            # makes sure required attributes are written into parameter file
            if MIN_LENGTH_ATTR not in config_sv:
                raise Exception("Min length must be specified on all SVs!")
            elif MAX_LENGTH_ATTR not in config_sv:
                raise Exception("Max length must be specified on all SVs!")
            elif TYPE_ATTR not in config_sv:
                raise Exception("\"Type\" attribute must be specified! For custom transformations, enter in \"Custom\"")

            # makes sure attributes are the correct datatype
            elif TYPE_ATTR in config_sv and not isinstance(config_sv[TYPE_ATTR], str):
                raise Exception("Invalid {} type for SV \'type\' attribute, str expected".format(type(config_sv[TYPE_ATTR])))

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
        try:
            config_entries = yaml.full_load(open(filename))
        except:
            raise Exception("YAML File {} failed to be open".format(filename))

        config = Config(**config_entries)

        return config

    def find_lcs(self, source_list, target_list, input_is_list=True):
        # Finds longest common subsequence between source_list and target_list such that if there is a dispersion event in source_list and target_list, it MUST belong in the lcs
        # If there is a dispersion (_) in source_list and target_list, then simply find the lcs of before and after the event and put them together
        # The lcs MUST have all the _ in source_list and target_list as a dispersion event cannot be transformed

        def lcs_list(list1, list2):
            # -> returns the longest common subsequence
            m, n = len(list1), len(list2)
        
            l = [[[] for i in range(n + 1)]
                for i in range(m + 1)]
        
            # bottom-up approach
            for i in range(m + 1):
                for j in range(n + 1):
                    if (i == 0 or j == 0):
                        l[i][j] = []
                    elif(list1[i - 1].upper() == list2[j - 1].upper()):
                        l[i][j] = l[i - 1][j - 1] + list(list2[j-1])
                    else:
                        l[i][j] = max(l[i - 1][j],
                                    l[i][j - 1], key=lambda lcs_entry: len(lcs_entry))
        
            return l[m][n]
        
        if not input_is_list:
            source_list = list(source_list)
            target_list = list(target_list)
        
        # remove placeholder from source_list and target_list
        source_list_new = [char for char in source_list if char != Symbols.PLACEHOLDER]
        target_list_new = [char for char in target_list if char != Symbols.PLACEHOLDER]

        # collects locations of dispersions (space) in source_list_new and target_list_new
        loc1 = [-1]
        loc1.extend([index for index in range(len(source_list_new)) if source_list_new[index].startswith(Symbols.DIS)])
        loc1.append(len(source_list_new))

        loc2 = [-1]
        loc2.extend([index for index in range(len(target_list_new)) if target_list_new[index].startswith(Symbols.DIS)])
        loc2.append(len(target_list_new))

        if len(loc1) != len(loc2):
            raise Exception("Unequal number of dispersion events ({}) across source_list_new and target_list_new, source_list_new: {}, target_list_new: {}".format(Symbols.DIS, source_list_new, target_list_new))
        
        common = list()
        for idx in range(1, len(loc1)):
            common.extend(lcs_list(source_list_new[loc1[idx-1] + 1: loc1[idx]], target_list_new[loc2[idx-1] + 1: loc2[idx]]))
            if idx != len(loc1) - 1:
                common.append(source_list_new[loc1[idx]])
            
        return common
    
    def export_to_bedpe(self, svs, bedfile, reset_file = True):
        def write_to_file(sv, source_chr, target_chr, source_s, source_e, target_s, target_e, transform, symbol, order = 0):
            assert (not symbol.startswith(Symbols.DIS))

            # do not write transformations of size 0
            if source_e > source_s or source_e == -1:
                with open(bedfile, "a") as fout:
                    row = [str(source_chr),
                            str(source_s),
                            str(source_e + 1),
                            str(target_chr),
                            str(target_s),
                            str(target_e + 1),
                            transform,
                            str(source_e - source_s),
                            str(int(sv.ishomozygous.value)) + "/1",
                            sv.name,
                            str(self.bedpe_counter),
                            str(order)]

                    fout.write("\t".join(row) + "\n")

        if reset_file:
            self.reset_file(bedfile)
        
        for x, sv in enumerate(svs):
            source, target = sv.source_unique_char, sv.target_unique_char
            events_dict = sv.events_dict

            # longest common subsequence represents the symbols that stay in place
            # the goal here is to minimize the number of transformations to get from source to target
            #source, target = "".join(source), "".join(target)
            same_place = self.find_lcs(source, target)
            #print(source, target, same_place)

            # Remove symbols in source which do not remain in same place
            index = 0
            for sr_symbol in source:
                if index >= len(same_place) or sr_symbol.upper() != same_place[index].upper():
                    event = events_dict[sr_symbol.upper()]
                    start, end = event.start, event.end
                    write_to_file(sv, event.source_chr, event.source_chr, start, end, start, start, Operations.DEL, sr_symbol)
                else:
                    index += 1
            
            # Add symbols from target not in same_place
            index = 0  # refers to position on same_place
            target_chr_index = 0
            target_chr = sv.target_block_chrs
            insert_point = sv.start
            order = 0
            for tr_symbol in target:
                # updates target chromosome to be exported
                # note this means that we've entered a new block
                if tr_symbol.startswith(Symbols.DIS):
                    target_chr_index += 1

                if index >= len(same_place) or tr_symbol.upper() != same_place[index].upper():   # current symbol does not remain in place
                    order += 1
                    event = events_dict[tr_symbol.upper()]
                    if tr_symbol.upper() in source:  # duplication
                        start, end = event.start, event.end
                        event_name = Operations.DUP
                    else:  # insertion detected
                        start, end = -1, -1
                        event_name = Operations.INS
                    write_to_file(sv, event.source_chr, target_chr[target_chr_index], start, end, insert_point, insert_point, event_name, tr_symbol, order)
                    
                    # Inversion applies in either case
                    if tr_symbol.islower():
                        write_to_file(sv, event.source_chr, target_chr[target_chr_index], start, end, insert_point, insert_point, Operations.INV, tr_symbol, order)
            
                elif tr_symbol.upper() == same_place[index].upper():   # current symbol does remain in place - the only possible event could be an inversion
                    event = events_dict[tr_symbol.upper()]
                    insert_point = event.end
                    order = 0
                    index += 1
                
                    # Inversion applies in either case
                    # Note that all dispersions should enter this case, but none will be exported because islower() returns False
                    if tr_symbol.islower():
                        start, end = event.start, event.end
                        write_to_file(sv, event.source_chr, target_chr[target_chr_index], start, end, start, end, Operations.INV, tr_symbol, order)
                
                else:
                    raise Exception("Target symbol {} does not fall into either case".format(tr_symbol))
            
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
    
    def close(self):
        # close all file handlers
        self.fin_export1.close()
        self.fin_export2.close()

class ErrorDetection():
    def __init__(self):
        pass

    @staticmethod
    def fail_if_any_overlapping(arr):
        # will raise Exception if any overlap between intervals is found
        # arr: list of tuples
        def is_overlapping(event_ranges, addition):
            # addition: tuple (start, end)
            # event_ranges: list containing tuples
            # checks if addition overlaps with any of the events already stored
    
            for event in event_ranges:
                if event[1] > addition[0] and event[0] < addition[1]:
                    return True
            return False
        for x, ele in enumerate(arr):
            if is_overlapping(arr[:x], ele):
                raise Exception("Overlapping Detected: {}".format(arr))
    
    @staticmethod
    def validate_symbols(transform):
        '''
        Ensures that transform has unique symbols other than dispersion events - specifically used for source sequence

        transform: str, source sequence
        '''
        present = dict()
        for symbol in transform:
            if symbol in present:
                raise Exception("Source transformation {} does not have unique symbols!".format(transform))
            elif symbol != Symbols.DIS:   # exclude dispersion events because they will always appear the same for user inputs
                present[symbol] = True
            

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
    
    '''
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



