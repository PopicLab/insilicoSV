import yaml
import sys
from constants import *
from pysam import FastaFile
import argparse
import os
import time

class Config():
    def __init__(self,**entries):
        self.__dict__.update(entries)
        self.run_checks()
        if "fail_if_placement_issues" not in self.__dict__:
            self.fail_if_placement_issues = False

        for config_sv in self.SVs:
            # handles cases where user enters length range for all components within SV or specifies different ranges
            if isinstance(config_sv["min_length"], int):
                config_sv["length_ranges"] = [(config_sv["min_length"], config_sv["max_length"])]
            else:
                config_sv["length_ranges"] = list(zip(config_sv["min_length"], config_sv["max_length"]))
            # make sure max_length >= min_length >= 0
            assert all(max_len >= min_len and min_len >= 0 for (min_len, max_len) in config_sv["length_ranges"]), "Max length must be >= min length for all SVs! Also ensure that all length values are >= 0."
            
            # supply filler values in for source and target if the type is custom
            config_sv["type"] = Variant_Type(config_sv["type"])
            if config_sv["type"] != Variant_Type.Custom:
                config_sv["source"] = None
                config_sv["target"] = None

    def run_checks(self):
        config_svs = self.SVs
        for config_sv in config_svs:
            # makes sure required attributes are written into parameter file
            if "min_length" not in config_sv:
                raise Exception("Min length must be specified on all SVs!")
            elif "max_length" not in config_sv:
                raise Exception("Max length must be specified on all SVs!")
            elif "type" not in config_sv:
                raise Exception("\"Type\" attribute must be specified! For custom transformations, enter in \"Custom\"")
            elif "number" not in config_sv:
                raise Exception("Number is a required parameter for all SVs")

            # makes sure attributes are the correct datatype
            elif "type" in config_sv and not isinstance(config_sv["type"], str):
                raise Exception("Invalid {} type for SV \'type\' attribute, str expected".format(type(config_sv["type"])))
            
class FormatterIO():
    def __init__(self, par_file):
        self.bedpe_counter = 1
        self.par_file = par_file

    def yaml_to_var_list(self):
        try:
            config_entries = yaml.full_load(open(self.par_file))
        except:
            raise Exception("YAML File {} failed to be open".format(self.par_file))

        config = Config(**config_entries)
        return config

    def export_to_bedpe(self, svs, bedfile, reset_file = True):
        '''
        Exports SVs changes to bedfile, uses target "blocks" to identify which subtype event is (INS, DUP, TRA, etc.)
        '''
        def write_to_file(sv, source_chr, source_s, source_e, target_chr, target_s, target_e, transform, symbol, event, order = 0):
            assert (not event.symbol.startswith(Symbols.DIS.value))

            # consider insertions
            if transform == Operations.INS.value:
                transform_length = event.length
                #transform_length = len(sv.events_dict[symbol].source_frag)
            else:
                transform_length = source_e - source_s

            # do not write events of size 0
            if event.length > 0:
                with open(bedfile, "a") as fout:
                    row = [str(source_chr),
                            str(source_s),
                            str(source_e + 1),
                            str(target_chr),
                            str(target_s),
                            str(target_e + 1),
                            transform,
                            str(transform_length),
                            str(int(sv.ishomozygous.value)) + "/1",
                            sv.name,
                            str(self.bedpe_counter),
                            str(order)]

                    fout.write("\t".join(row) + "\n")
        def symbol_is_inversion(symbol):
            return any(c.islower() for c in symbol)

        if reset_file:
            self.reset_file(bedfile)
        
        for x, sv in enumerate(svs):
            source, target = sv.source_unique_char, sv.target_unique_char
            encoding = sv.events_dict    # maps symbol like A or B to base pairs on reference 
            #print("Encode_dict: ", encoding)

            assert (sv.start != None and sv.end != None) # start & end should have been defined alongside event positions
            curr_pos = sv.start
            curr_chr = sv.start_chr

            for idx, block in enumerate(sv.target_symbol_blocks):
                order = 0
                for x, symbol in enumerate(block):
                    event = encoding[symbol.upper()[0]]

                    # duplication/inverted-duplications
                    if len(symbol) > 1 and symbol[1] == Symbols.DUP_MARKING.value:
                        if symbol_is_inversion(symbol):
                            event_name = Operations.INVDUP.value
                        else:
                            event_name = Operations.DUP.value
                        # consider dispersed duplications
                        if event.end != curr_pos or event.source_chr != curr_chr:
                            event_name = "d" + event_name
                        
                        order += 1
                        write_to_file(sv, event.source_chr, event.start, event.end, curr_chr, curr_pos, curr_pos, event_name, event.symbol, event, order)
                    
                    else:  # symbol is an original
                        # insertions
                        if symbol[0].upper() not in sv.source_unique_char:
                            order += 1
                            write_to_file(sv, curr_chr, curr_pos, curr_pos, curr_chr, curr_pos, curr_pos, Operations.INS.value, event.symbol, event, order)
                        
                        # translocations - original symbol
                        elif len(symbol) == 1 and idx != event.original_block_idx:
                            order += 1
                            if symbol_is_inversion(symbol):
                                event_name = Operations.INVTRA.value
                            else:
                                event_name = Operations.TRA.value
                            write_to_file(sv, event.source_chr, event.start, event.end, curr_chr, curr_pos, curr_pos, event_name, event.symbol, event, order)

                        # inversions
                        elif symbol_is_inversion(symbol):
                            order = 0
                            curr_pos = event.end       # this symbol was already in source
                            write_to_file(sv, event.source_chr, event.start, event.end, curr_chr, event.start, event.end, Operations.INV.value, event.symbol, event, order)
                        
                        # identity - original symbol did not change or move
                        else:
                            order = 0
                            curr_pos = event.end

                # find dispersion event right after block to update position and chromosome to edit on
                if idx < len(sv.target_symbol_blocks) - 1:
                    dis_event = encoding[Symbols.DIS.value + str(idx + 1)]  # find the nth dispersion event
                    curr_pos = dis_event.end
                    curr_chr = dis_event.source_chr
            
            # deletions - any original symbols not detected in target sequence are deleted
            for symbol in source:
                # do not delete dispersion events or symbols already in target
                if not symbol.startswith(Symbols.DIS.value) and symbol not in target and symbol.lower() not in target:
                    event = encoding[symbol]
                    order = 0
                    write_to_file(sv, event.source_chr, event.start, event.end, event.source_chr, event.start, event.start, Operations.DEL.value, event.symbol, event, order)
            self.bedpe_counter += 1


    def export_piece(self, id, edits, fasta_out, fasta_file, verbose = False):
        '''
        Exports list of changes from simulator to fasta file

        id: chr_id to apply edits to
        edits: list with interval and new_frag, replace reference at the interval with new_frag
        fasta_out: Fasta file to export changes to
        fasta_file: FastaFile with access to reference
        '''    
        def get_buffer_size(curr_pos, target):
            # returns the number of bases to read at the given position
            # max bases to take in is MAX_BUFFER_SIZE
            if target - curr_pos > MAX_BUFFER_SIZE:
                return MAX_BUFFER_SIZE
            else:
                return target - curr_pos

        with open(fasta_out,"a") as fout_export:
            if id not in fasta_file.references:
                raise KeyError("ID {} not found in inputted fasta file".format(id))
            
            if verbose:
                print("New ID: ", id)

            # Write in chr_id
            fout_export.write(">" + str(id) + "\n")
            chr_variants = list(edits)   # given as (start, end, new_frag)
            chr_variants.sort()
            chr_variants.append([fasta_file.get_reference_length(id), fasta_file.get_reference_length(id), ""]) # to ensure that all bases after last variant are included

            # apply changes by traversing reference from start to end - this is possible since chr_variants is sorted
            pos = 0
            for variant in chr_variants:
                var_start, var_end = variant[0], variant[1]
                # write all reference bases up until interval, which indicates that a change will occur
                while pos < var_start:
                    appropriate_buffer = get_buffer_size(pos, var_start)
                    c = fasta_file.fetch(id, pos, pos + appropriate_buffer)
                    fout_export.write(c)
                    pos += appropriate_buffer

                # add in replacement fragment
                assert (pos == var_start), "Replacement fragment about to be inserted at position {} instead of var_start {}".format(pos, var_start)
                # after var_start, loop will ignore up till var_end
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

def collect_args():
    parser = argparse.ArgumentParser(description='insilicoSV is a software to design and simulate complex structural variants, both novel and known.')
    parser.add_argument("ref", help="FASTA reference file")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("prefix", help="Used for naming haplotype, bed, and stats file")
    #parser.add_argument("hap1", help="First output FASTA file (first haplotype)")
    #parser.add_argument("hap2", help = "Second output FASTA file (second haplotype)")
    #parser.add_argument("bedpe", help = "BEDPE file location to store variant info")
    #parser.add_argument("stats", help = "File location for stats file")
    parser.add_argument("-r", "--root", action="store", metavar="DIR", dest="root_dir", help="root directory for all files given")

    args = parser.parse_args()
    #io = [args.ref, args.config, args.hap1, args.hap2, args.bedpe, args.stats]
    args_dict = {"ref": args.ref, "config": args.config, "hap1": args.prefix + ".hapA.fa",
                "hap2": args.prefix + ".hapB.fa", "bedpe": args.prefix + ".bed", "stats": args.prefix + ".stats.txt"}
    if args.root_dir:
        for key, curr_path in args_dict.items():
            args_dict[key] = os.path.join(args.root_dir, curr_path)
    return args_dict



