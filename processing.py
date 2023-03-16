import yaml
import sys

import constants
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
        # optional config feature storing the path to a overlap_events bed file
        if "overlap_events" in entries:
            self.__dict__["overlap_events"] = entries["overlap_events"]
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
        valid_keys = ["sim_settings", "SVs", "overlap_events"]  # valid arguments at the top level
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

    def write_to_file(self, sv, bedfile, source_s, source_e, target_s, target_e, transform, event, nth_sv, order=0):
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
                       str(source_e),
                       str(event.source_chr),
                       str(target_s),
                       str(target_e),
                       transform,
                       str(transform_length),
                       str(int(sv.ishomozygous.value)) + "/1",
                       sv.name,
                       str(nth_sv),
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
    def get_event_target_operation(ev, target_events_dict, source_events_dict):
        # helper function to return target_intvl and operation for multi-source events
        # need to enumerate the possible modifications to set the right operation
        # A -> A'
        if ev + Symbols.DUP_MARKING.value in target_events_dict.keys():
            trg_sym = ev + Symbols.DUP_MARKING.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), \
                   Operations.DUP.value if ev in target_events_dict.keys() else Operations.TRA.value
        # A -> a'
        elif ev.lower() + Symbols.DUP_MARKING.value in target_events_dict.keys():
            trg_sym = ev.lower() + Symbols.DUP_MARKING.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INVDUP.value
        # A -> a
        elif ev.lower() in target_events_dict.keys():
            trg_sym = ev.lower()
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INV.value
        # A -> A* (in the case of a custom event in which an event is divergently duplicated)
        elif ev + Symbols.DIV_MARKING.value in target_events_dict.keys():
            trg_sym = ev + Symbols.DIV_MARKING.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.DIV.value
        # A -> A (insertion if source A is undefined, identity otherwise)
        elif ev in target_events_dict.keys():
            return (target_events_dict[ev].start, target_events_dict[ev].end), \
                Operations.INS.value if source_events_dict[ev].start is None else Operations.IDENTITY.value
        # A -> [none]
        elif ev not in [sym[0] for sym in target_events_dict.keys()]:
            return (source_events_dict[ev].start, source_events_dict[ev].end), Operations.DEL.value
        # otherwise unknown mapping
        else:
            return (source_events_dict[ev].start, source_events_dict[ev].end), Operations.UNDEFINED.value

    @staticmethod
    def postprocess_record_params(sv, sv_record_info):
        # arrange the bed_record parameter dictionaries in order of ascending source interval start position
        # and assign order values to the relevant entries
        # --> for TRA/INS/DUP events with the same target position, 'order' describes the order in which they
        # --> are compiled (i.e., the order in which they appear in the target sequence)
        order = 0
        ins_pos = None
        for block in sv.target_symbol_blocks:
            for target_event in block:
                if target_event.symbol.startswith(Symbols.DIS_MARKING.value) or \
                        target_event.symbol in sv_record_info.keys():  # <- prevent collision with A' and A if both in target
                    continue
                src_sym = target_event.symbol[0].upper()
                if sv_record_info[src_sym]['transform'] in constants.NONZERO_ORDER_OPERATIONS:
                    if ins_pos is None:
                        ins_pos = sv_record_info[src_sym]['target_s']
                        order += 1
                    elif sv_record_info[src_sym]['target_s'] == ins_pos:
                        order += 1
                else:
                    ins_pos = None
                    order = 0
                sv_record_info[src_sym]['order'] = order
        return sorted([params for params in sv_record_info.values()], key=lambda params: params['source_s'])

    def export_to_bedpe(self, svs, bedfile, ins_fasta=None, reset_file=True):
        if reset_file:
            utils.reset_file(bedfile)
            if ins_fasta:
                utils.reset_file(ins_fasta)
        nth_sv = 0
        for sv in svs:
            nth_sv += 1
            # SVs with multiple source events will be split into multiple bed records (one for each)
            if len(sv.events_dict) == 1:
                ev = list(sv.sv_blocks.target_events_dict.values())[0] if sv.type == Variant_Type.INS\
                        else list(sv.events_dict.values())[0]
                op = self.get_event_target_operation(ev.symbol, sv.sv_blocks.target_events_dict, sv.events_dict)[1]
                record_info = {'source_s': ev.start, 'source_e': ev.end, 'target_s': ev.start, 'target_e': ev.end,
                               'transform': op, 'sv': sv, 'event': ev, 'bedfile': bedfile, 'nth_sv': nth_sv,
                               # for simple events, order cant be > 1 and only depends on event type
                               'order': int(op in constants.NONZERO_ORDER_OPERATIONS)}
                self.write_to_file(**record_info)
                # simple INS -> option to write novel inserted sequences to separate .fa at ins_fasta
                if op == Operations.INS.value:
                    self.export_insertions(sv.start_chr, ev.start, ev.source_frag, ins_fasta)
            else:
                # multiple source events: source intervals taken from the source events
                # and target intervals taken from corresponding target events (if no match, then deletion)
                sv_record_info = {}
                for ev in sv.events_dict.values():
                    # skip dispersion events
                    if ev.symbol.startswith(Symbols.DIS_MARKING.value):
                        continue
                    # --> target intvl and transform will be determined in each of the specific cases of how the
                    # --> source event is mapped into the target
                    sv_record_info[ev.symbol] = {'source_s': ev.start, 'source_e': ev.end, 'sv': sv, 'event': ev, 'bedfile': bedfile, 'nth_sv': nth_sv}
                    (target_s, target_e), operation = self.get_event_target_operation(ev.symbol, sv.sv_blocks.target_events_dict, sv.events_dict)
                    sv_record_info[ev.symbol]['target_s'] = target_s
                    sv_record_info[ev.symbol]['target_e'] = target_e
                    sv_record_info[ev.symbol]['transform'] = operation
                for param_dict in self.postprocess_record_params(sv, sv_record_info):
                    self.write_to_file(**param_dict)

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
        vcf_file.header.info.add('OVERLAP_EV', number=1, type='String', description="Bool. indicator for the event being"
                                                                                    "placed at an overlap_events interval")
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
            if sv.overlap_event is not None:
                info_field['OVERLAP_EV'] = 'True'

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


