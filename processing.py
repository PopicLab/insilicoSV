import yaml
import sys
from constants import *
from pysam import FastaFile
import pysam
import argparse
import os
import time
import utils

class Config():
    def __init__(self,**entries):
        self.__dict__.update(DEFAULT_CONFIG)
        self.__dict__.update({"SVs": entries["SVs"]})
        if "sim_settings" in entries:   # if the user changed any of the default settings for the sim, they get updated here
            self.__dict__["sim_settings"].update(entries["sim_settings"])
        # optional config feature storing the path to a repeatmasker bed file
        if "repeatmasker" in entries:
            self.__dict__["repeatmasker"] = entries["repeatmasker"]
        # artificial "keys" attribute to access the keys of config.__dict__
        self.keys = self.__dict__.keys()

        if "vcf_path" not in self.__dict__["SVs"][0]:
            self.run_checks_randomized(entries)
            for config_sv in self.SVs:
                if "avoid_intervals" in config_sv:
                    continue
                # handles cases where user enters length range for all components within SV or specifies different ranges
                if isinstance(config_sv["min_length"], int):
                    config_sv["length_ranges"] = [(config_sv["min_length"], config_sv["max_length"])]
                else:
                    config_sv["length_ranges"] = list(zip(config_sv["min_length"], config_sv["max_length"]))
                # make sure max_length >= min_length >= 0
                assert all(max_len >= min_len >= 0 for (min_len, max_len) in config_sv["length_ranges"]), "Max length must be >= min length for all SVs! Also ensure that all length values are >= 0."

                # use Enum for variant type and nonvariant event type
                if Nonvariant_Event_Type.has_value(config_sv["type"]):
                    config_sv["type"] = Nonvariant_Event_Type(config_sv["type"])
                else:
                    config_sv["type"] = Variant_Type(config_sv["type"])
                if config_sv["type"] != Variant_Type.Custom:
                    config_sv["source"] = None
                    config_sv["target"] = None

    def run_checks_randomized(self, entries):
        '''
        check method for yaml given with SVs given for randomized placement on reference
        '''
        # entries: dict representing full config file parsed from yaml
        config_svs = self.SVs
        for config_sv in config_svs:
            # "avoid_intervals" the optional argument for specifying a vcf whose event intervals
            # will be avoided during random placement of simulated events
            if "avoid_intervals" in config_sv:
                continue
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
        valid_optional_par = ["fail_if_placement_issues", "max_tries", "generate_log_file", "filter_small_chr", "prioritize_top"] # valid arguments within sim_settings
        for parameter in self.sim_settings:
            if parameter not in valid_optional_par:
                raise Exception("\"{}\" is an invalid argument under sim_settings".format(parameter))
        valid_keys = ["sim_settings", "SVs", "repeatmasker"]  # valid arguments at the top level
        for key in entries:
            if key not in valid_keys:
                raise Exception("Unknown argument \"{}\"".format(key))

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

    def write_to_file(self, sv, bedfile, source_s, source_e, target_s, target_e, transform, event, order=0):
        assert (not event.symbol.startswith(Symbols.DIS_MARKING.value))
        # consider insertions
        if transform == Operations.INS.value:
            transform_length = event.length
        else:
            transform_length = source_e - source_s

        # do not write events of size 0
        if event.length > 0:
            with open(bedfile, "a") as fout:
                row = [str(event.source_chr),
                       str(source_s),
                       str(source_e + 1),
                       str(event.source_chr),
                       str(target_s),
                       str(target_e + 1),
                       transform,
                       str(transform_length),
                       str(int(sv.ishomozygous.value)) + "/1",
                       sv.name,
                       str(self.bedpe_counter),
                       str(order)]

                fout.write("\t".join(row) + "\n")

    @staticmethod
    def symbol_is_inversion(symbol):
        return any(c.islower() for c in symbol)

    @staticmethod
    def export_insertions(chr, start_pos, seq, ins_fasta):
        '''
        Exports foreign insertion sequences to separate fasta file, append only
        '''
        with open(ins_fasta, "a") as fout_ins:
            fout_ins.write(">{}_{}\n".format(chr, start_pos))
            fout_ins.write("{}\n".format(seq))

    @staticmethod
    def get_composite_sv_event_info(ev, target_events_dict, source_events_dict):
        # helper function to return target_intvl and operation for multi-source events
        # need to enumerate the possible modifications to set the right operation
        # A -> A'
        if ev + Symbols.DUP_MARKING.value in target_events_dict.keys():
            trg_sym = ev + Symbols.DUP_MARKING.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.DUP.value
        # A -> a'
        elif ev.lower() + Symbols.DUP_MARKING.value in target_events_dict.keys():
            trg_sym = ev.lower() + Symbols.DUP_MARKING.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INVDUP.value
        # A -> a
        elif ev.lower() in target_events_dict.keys():
            trg_sym = ev.lower()
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INV.value
        # A -> A* (in the case of a custom event in which an event is divergently duplicated)
        elif ev + Symbols.DIV_MARKING in target_events_dict.keys():
            trg_sym = ev + Symbols.DIV_MARKING
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.DIV.value
        # otherwise unknown mapping
        else:
            return (source_events_dict[ev].start, source_events_dict[ev].end), Operations.UNDEFINED.value

    def export_to_bedpe(self, svs, bedfile, ins_fasta, reset_file=True):
        if reset_file:
            utils.reset_file(bedfile)
            utils.reset_file(ins_fasta)
        ## for reference, the input arguments to write_to_file:
        # sv, bedfile, source_s, source_e, target_s, target_e, transform, symbol, event, order = 0
        for sv in svs:
            # create target events dict for lookup of corresponding source/target events within SV
            target_events_dict = sv.sv_blocks.target_events_dict()
            # SVs with multiple source events will be split into multiple bed records (one for each)
            if len(sv.events_dict) < 2:
                # simple DEL
                if len(sv.events_dict) == 1:
                    ev = list(sv.events_dict.values())[0]
                    op = Operations.DEL.value
                # simple INS
                else:
                    ev = list(target_events_dict.values())[0]
                    op = Operations.INS.value
                record_info = {'source_s': ev.start,
                               'source_e': ev.end,
                               'target_s': ev.start,
                               'target_e': ev.end,
                               'transform': op,
                               'sv': sv,
                               'event': ev,
                               'bedfile': bedfile}
                self.write_to_file(**record_info)
            else:
                # TODO: add logic to calculate 'order' value for each relevant composite event
                # multiple source events: source intervals taken from the source events
                # and target intervals taken from corresponding target events (if no match, then deletion)
                # --> dict keyed on source symbol with values giving source and target intervals
                sv_record_info = {}
                for ev in sv.events_dict.keys():
                    # --> target intvl and transform will be determined in each of the specific cases of how the
                    # --> source event is mapped into the target
                    sv_record_info[ev] = {'source_s': sv.events_dict[ev].start,
                                          'source_e': sv.events_dict[ev].end,
                                          'sv': sv,
                                          'event': ev,
                                          'bedfile': bedfile}
                    # source event appears unaltered in target: identity subevent
                    if ev in target_events_dict.keys():
                        sv_record_info[ev]['target_s'] = target_events_dict[ev].start
                        sv_record_info[ev]['target_e'] = target_events_dict[ev].end
                        sv_record_info[ev]['transform'] = Operations.IDENTITY.value
                    # source event appears altered in target: dup, inv, etc.
                    elif ev in [sym[0].upper() for sym in target_events_dict.keys()]:
                        (target_s, target_e), operation = self.get_composite_sv_event_info(ev, target_events_dict, sv.events_dict)
                        sv_record_info[ev]['target_s'] = target_s
                        sv_record_info[ev]['target_e'] = target_e
                        sv_record_info[ev]['transform'] = operation
                    # source event doesn't appear in any form in target events (deletion)
                    else:
                        sv_record_info[ev]['target_s'] = sv.events_dict[ev].start
                        sv_record_info[ev]['target_e'] = sv.events_dict[ev].end
                        sv_record_info[ev]['transform'] = Operations.DEL.value
                # --> write bed records from interval_dict contents in the order of source interval start
                for param_dict in sorted([params for params in sv_record_info.values()], key=lambda params: params.source_s):
                    self.write_to_file(**param_dict)

    def export_to_bedpe_OLD(self, svs, bedfile, ins_fasta, reset_file = True):
        '''
        Exports SVs changes to bedfile and foreign insertions to separate fasta file, uses target "blocks" to identify which subtype event is (INS, DUP, TRA, etc.)
        svs: list of Structural Variant objects
        bedfile: File location for bedpe file
        ins_fasta: File location to export foreign insertions
        reset_file: True = delete any previous content in bed file, False = append
        '''

        if reset_file:
            utils.reset_file(bedfile)
            utils.reset_file(ins_fasta)

        for x, sv in enumerate(svs):
            source, target = sv.source_unique_char, sv.target_unique_char
            encoding = sv.events_dict    # maps symbol like A or B to base pairs on reference

            assert (sv.start != None and sv.end != None) # start & end should have been defined alongside event positions
            curr_pos = sv.start
            curr_chr = sv.start_chr
            for idx, block in enumerate(sv.target_symbol_blocks):
                order = 0
                for x, ev in enumerate(block):
                    symbol = ev.symbol
                    event = encoding[symbol.upper()[0]]

                    # duplication/inverted-duplications
                    if len(symbol) > 1 and symbol[1] == Symbols.DUP_MARKING.value:
                        if self.symbol_is_inversion(symbol):
                            event_name = Operations.INVDUP.value
                        else:
                            event_name = Operations.DUP.value
                        # consider dispersed duplications
                        if event.end != curr_pos or event.source_chr != curr_chr:
                            event_name = "d" + event_name
                        order += 1
                        self.write_to_file(sv, bedfile, event.source_chr, event.start, event.end, curr_chr, curr_pos, curr_pos, event_name, event.symbol, event, order)

                    else:  # symbol is an original
                        # insertions
                        if symbol[0].upper() not in sv.source_unique_char:
                            order += 1
                            self.write_to_file(sv, bedfile, curr_chr, curr_pos, curr_pos, curr_chr, curr_pos, curr_pos, Operations.INS.value, event.symbol, event, order)
                            self.export_insertions(curr_chr, curr_pos, event.source_frag, ins_fasta)   # foreign insertion sequences need to be exported separate from bed file

                        # translocations - original symbol
                        # elif len(symbol) == 1 and idx != event.original_block_idx:
                        #     order += 1
                        #     if self.symbol_is_inversion(symbol):
                        #         event_name = Operations.INVTRA.value
                        #     else:
                        #         event_name = Operations.TRA.value
                        #     self.write_to_file(sv, bedfile, event.source_chr, event.start, event.end, curr_chr, curr_pos, curr_pos, event_name, event.symbol, event, order)

                        # inversions
                        elif self.symbol_is_inversion(symbol):
                            order = 0
                            curr_pos = event.end       # this symbol was already in source
                            self.write_to_file(sv, bedfile, event.source_chr, event.start, event.end, curr_chr, event.start, event.end, Operations.INV.value, event.symbol, event, order)

                        # identity - original symbol did not change or move
                        else:
                            order = 0
                            curr_pos = event.end

                # find dispersion event right after block to update position and chromosome to edit on
                if idx < len(sv.target_symbol_blocks) - 1:
                    dis_event = encoding[Symbols.DIS_MARKING.value + str(idx + 1)]  # find the nth dispersion event
                    curr_pos = dis_event.end
                    curr_chr = dis_event.source_chr

            # deletions - any original symbols not detected in target sequence are deleted
            for symbol in source:
                # do not delete dispersion events or symbols already in target
                if not symbol.startswith(Symbols.DIS_MARKING.value) and symbol not in target and symbol.lower() not in target:
                    event = encoding[symbol]
                    order = 0
                    self.write_to_file(sv, bedfile, event.source_chr, event.start, event.end, event.source_chr, event.start, event.start, Operations.DEL.value, event.symbol, event, order)
            self.bedpe_counter += 1

    def export_to_vcf(self, svs, stats, vcffile):
        with open(vcffile, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            for chrm, chrm_len in stats.chr_lengths.items():
                vcf.write("##contig=<ID=%s,length=%d>\n" % (chrm, chrm_len))
            vcf.write("#%s\n" % "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                                           "SAMPLE"]))
        # *** This will throw an error with pysam version 0.18, need 0.16.0.1
        vcf_file = pysam.VariantFile(vcffile)
        vcf_file.header.info.add('END', number=1, type='Integer', description="End position of the variant "
                                                                              "described in this record")
        vcf_file.header.info.add('CIPOS', number=2, type='Integer', description="Confidence interval around POS for "
                                                                                "imprecise variants")
        vcf_file.header.info.add('CIEND', number=2, type='Integer', description="Confidence interval around END for "
                                                                                "imprecise variants")
        vcf_file.header.info.add('SVTYPE', number=1, type='String', description="Type of structural variant")
        vcf_file.header.info.add('SVLEN', number=1, type='Integer', description="Length of structural variant")
        vcf_file.header.info.add('SVMETHOD', number=1, type='String', description="SV detection method")
        vcf_file.header.info.add('TARGET', number=1, type='Integer', description="Target location for divergent repeat")
        vcf_file.header.info.add('DIV_REPEAT', number=1, type='String', description="Divergent repeat segment places at target locus")
        vcf_file.header.info.add('RM_OVERLAP', number=1, type='String', description="Bool. indicator for the event being"
                                                                                    "placed at a repeatmasker interval")
        vcf_file.header.formats.add('GT', number=1, type='String', description="Genotype")

        vcf_out_file = pysam.VariantFile(vcffile, 'w', header=vcf_file.header)

        for sv in svs:
            zyg = (1, 1) if sv.ishomozygous == Zygosity.HOMOZYGOUS else (0, 1)
            dispersion_target = None
            if sv.type in DISPERSION_TYPES:
                # --> going to read the source and target intervals off the start/end positions of the events dict
                source_event = sv.events_dict['A']
                disp_event = sv.events_dict['_1']
                rec_start = source_event.start
                rec_end = source_event.end
                # setting TARGET to the correct end of the dispersion
                if disp_event.start == source_event.end:
                    dispersion_target = disp_event.end
                else:
                    dispersion_target = disp_event.start
            else:
                # start/end given by the min/max changed fragment interval positions
                rec_start = min([frag[1] for frag in sv.changed_fragments])
                rec_end = max(frag[2] for frag in sv.changed_fragments)
            if dispersion_target is not None:
                info_field = {'SVTYPE': sv.type.value, 'SVLEN': rec_end - rec_start, 'TARGET': dispersion_target}
            else:
                if sv.type == Variant_Type.INS:
                    # special case of simple INS: sv length \neq (sv end - sv start)
                    info_field = {'SVTYPE': sv.type.value, 'SVLEN': sv.events_dict['A'].length}
                else:
                    info_field = {'SVTYPE': sv.type.value, 'SVLEN': rec_end - rec_start}
            if sv.repeatmasker_event is not None:
                info_field['RM_OVERLAP'] = 'True'

            vcf_record = vcf_out_file.header.new_record(contig=sv.start_chr, start=rec_start, stop=rec_end,
                                                        alleles=['N', sv.type.value], id=sv.type.value,
                                                        info=info_field,
                                                        qual=100, filter='PASS',
                                                        samples=[{'GT': zyg}])
            vcf_out_file.write(vcf_record)

        vcf_out_file.close()


    def export_variants_to_fasta(self, id, edits, fasta_out, fasta_file, verbose = False):
        '''
        Exports list of changes from simulator to fasta file

        id: chr_id to apply edits to
        edits: list with interval and new_frag, replace reference at the interval with new_frag
        fasta_out: Fasta file to export changes to
        fasta_file: FastaFile with access to reference
        '''
        with open(fasta_out,"a") as fout_export:
            if id not in fasta_file.references:
                raise KeyError("ID {} not found in inputted fasta file".format(id))

            if verbose:
                print("New ID: ", id)

            # Write in chr_id
            fout_export.write(">" + str(id) + "\n")
            chr_variants = list(edits)   # given as (start, end, new_frag)
            chr_variants.sort()
            # *** I think this is colliding with the case of a deletion of the entire input ref sequence
            chr_variants.append([fasta_file.get_reference_length(id), fasta_file.get_reference_length(id), ""]) # to ensure that all bases after last variant are included

            # apply changes by traversing reference from start to end - this is possible since chr_variants is sorted
            pos = 0
            for variant in chr_variants:
                var_start, var_end = variant[0], variant[1]
                # write all reference bases up until interval, which indicates that a change will occur
                while pos < var_start:
                    appropriate_buffer = MAX_BUFFER_SIZE if var_start - pos > MAX_BUFFER_SIZE else var_start - pos
                    c = fasta_file.fetch(id, pos, pos + appropriate_buffer)
                    fout_export.write(c)
                    pos += appropriate_buffer

                # add in replacement fragment
                assert (pos == var_start), "Replacement fragment about to be inserted at position {} instead of var_start {}".format(pos, var_start)
                fout_export.write(variant[2])

                pos = var_end
            fout_export.write("\n")

    def close(self):
        # close all file handlers
        self.fin_export1.close()
        self.fin_export2.close()

def collect_args():
    parser = argparse.ArgumentParser(description='insilicoSV is a software to design and simulate complex structural variants, both novel and known.')
    parser.add_argument("ref", help="FASTA reference file")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("prefix", help="Used for naming all output files (haplotype, bed, and stats)")
    parser.add_argument("-r", "--root", action="store", metavar="DIR", dest="root_dir", help="root directory for all files given")

    args = parser.parse_args()
    args_dict = {"ref": args.ref, "config": args.config, "ins_fasta": args.prefix + ".insertions.fa", "hap1": args.prefix + ".hapA.fa",
                "hap2": args.prefix + ".hapB.fa", "bedpe": args.prefix + ".bed", "stats": args.prefix + ".stats.txt", "log_file": args.prefix + ".log"}
    if args.root_dir:
        for key, curr_path in args_dict.items():
            args_dict[key] = os.path.join(args.root_dir, curr_path)
    return args_dict


