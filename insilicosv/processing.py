import argparse
import os
import sys
import time

from pysam import FastaFile
import pysam
import yaml

from insilicosv import utils
from insilicosv.constants import \
    MAX_BUFFER_SIZE, Variant_Type, DISPERSION_TYPES, Operations, \
    NONZERO_ORDER_OPERATIONS, Symbols, SV_KEY

class FormatterIO:
    def __init__(self, par_file):
        self.bedpe_counter = 1
        self.par_file = par_file
        self.config = None

    @staticmethod
    def run_checks_randomized(config):
        """
        check method for yaml given with SVs given for randomized placement on reference
        """
        config_svs = config['variant_sets']
        for config_sv in config_svs:
            if "blacklist_regions" in config_sv:
                continue
            elif "type" not in config_sv:
                raise Exception("\"Type\" attribute must be specified! For custom transformations, enter in \"Custom\"")
            elif config_sv["type"] == "SNP":  # SNP events are only specified by count (size is deterministic)
                if "number" in config_sv and isinstance(config_sv["number"], int) and config_sv["number"] > 0:
                    continue
                else:
                    raise Exception("Number (of type int > 0) is a required parameter for all SVs")
            if "length_ranges" not in config_sv:
                raise Exception("Length ranges must be specified on all SVs!")
            if "number" not in config_sv:
                raise Exception("Number is a required parameter for all SVs")

            elif "type" in config_sv and not isinstance(config_sv["type"], str):
                raise Exception("Invalid {} type for SV \'type\' attribute, str expected".format(type(config_sv["type"])))
        valid_optional_par = ["fail_if_placement_issues", "max_tries", "generate_log_file", "filter_small_chr",
                              "prioritize_top", "homozygous_only", "reference", "min_intersv_dist"]  # valid arguments within sim_settings
        for parameter in config['sim_settings']:
            if parameter not in valid_optional_par:
                raise Exception("\"{}\" is an invalid argument under sim_settings".format(parameter))
        valid_keys = ["sim_settings", "variant_sets", "overlap_regions", "blacklist_regions"]  # valid arguments at the top level
        for key in config:
            if key not in valid_keys:
                raise Exception("Unknown argument \"{}\"".format(key))

    def postproc_config_dict(self):
        if 'sim_settings' not in self.config:
            raise Exception("Must include \'sim_settings\' sections specifying at least \'reference\' path")
        if "filter_small_chr" in self.config and not isinstance(self.config["filter_small_chr"], int):
            raise Exception("Must provide value of type int to \'filter_small_chr\'")
        if "reference" not in self.config["sim_settings"]:
            raise Exception("Must include reference FASTA file in \'reference\' field of \'sim_settings\'")
        elif self.config["sim_settings"]["reference"].split(".")[-1] not in ["fa", "fna", "fasta"]:
            raise Exception("Input reference must be of type .fa, .fna, or .fasta")
        if "vcf_path" not in self.config["variant_sets"][0]:
            self.run_checks_randomized(self.config)
        for config_sv in self.config['variant_sets']:
            if "vcf_path" in config_sv:
                continue
            assert "type" in config_sv, "\'type\' missing from \'variant_sets\' config entry"
            if any(ovlp_field in config_sv for ovlp_field in ["overlap_type", "overlap_region_type", "overlap_component"]):
                assert "overlap_regions" in self.config, "Must provide \'overlap_regions\' entry if \'overlap_type\', " \
                                                       "\'overlap_region_type\', or \'overlap_component\' are specified"
            config_sv["type"] = Variant_Type(config_sv["type"])
            if config_sv["type"] != Variant_Type.Custom:
                config_sv["source"] = None
                config_sv["target"] = None
            if 'overlap_component' in config_sv and config_sv['overlap_component'] in ['source', 'target']:
                if not (config_sv['type'] in DISPERSION_TYPES or (
                        config_sv['source'] is not None and Symbols.DIS.value in config_sv['source'])):
                    raise Exception(
                        'Cannot specify \'overlap_component\': {\'source\'/\'target\'} for non-dispersion SVs')
            elif 'overlap_component' in config_sv and config_sv['overlap_component'] == 'full_sv':
                if config_sv['type'] in DISPERSION_TYPES or (
                        config_sv['source'] is not None and Symbols.DIS.value in config_sv['source']):
                    raise Exception('Cannot specify \'overlap_component\': \'full_sv\' for dispersion SVs')
            frags = self.get_grammar(config_sv)
            unbounded_disps = []
            # SV event length specification - not applicable for SNPs
            if config_sv["type"] != Variant_Type.SNP:
                if not isinstance(config_sv["length_ranges"], list):
                    raise Exception("Must provide list of tuples to \'length\'")
                if len(config_sv["length_ranges"]) != len(frags):
                    raise Exception(f"SV of type {config_sv['type']} requires {len(frags)} length_ranges")
                for i, frag in enumerate(frags):
                    min_len, max_len = config_sv["length_ranges"][i]
                    if frag[0] != '_':
                        assert min_len is not None and max_len is not None, "Min and max length must be non-None for non-dispersion fragments"
                        assert max_len >= min_len >= 0, "Max length must be >= min length for all SVs! Also ensure that all length values are >= 0."
                    else:
                        assert min_len is not None, "Min length must be non-None for dispersion fragments"
                        if max_len is not None:
                            assert max_len >= min_len >= 0, "Max length must be >= min length for all SVs! Also ensure that all length values are >= 0."
                        else:
                            unbounded_disps.append(i)
                            assert min_len >= 0, "Min length must be >= 0 for all SVs!"
            if "divergence_prob" in config_sv:
                assert isinstance(config_sv["divergence_prob"], (int, float)), \
                    "Must give \'divergence_prob\'"
                assert 1 >= config_sv["divergence_prob"] >= 0, "divergence_prob must be in [0,1]"
                assert config_sv["type"] != "SNP" or config_sv["divergence_prob"] == 1, \
                    "divergence_prob for SNPs is always 1"

            # if overlap_region_type given as list, must provide a list of equal length to the number of frags in SV
            if "overlap_region_type" in config_sv and isinstance(config_sv["overlap_region_type"], list):
                assert len(frags) == len(config_sv["overlap_region_type"]), "if \'overlap_region_type\' given as list, " \
                                                                           "must have length equal to number of SV components"
                # Can only support multiple non-None overlap frags if separated by unbounded dispersions
                roi_idxs = [i for i, roi in enumerate(config_sv['overlap_region_type']) if roi is not None]
                if len(roi_idxs) > 1:
                    for i in range(len(roi_idxs) - 1):
                        assert any([roi_idxs[i] < elt <= roi_idxs[i + 1] for elt in unbounded_disps]), \
                            '\'overlap_region_type\' can only be non-None for multiple fragments that are separated by an unbounded dispersion'

        # setting default values for sim_settings fields
        if 'max_tries' not in self.config['sim_settings']:
            self.config['sim_settings']['max_tries'] = 50
        if 'fail_if_placement_issues' not in self.config['sim_settings']:
            self.config['sim_settings']['fail_if_placement_issues'] = True

        # assign a unique ID to each variant set
        for variant_set_id, variant_set in enumerate(self.config['variant_sets']):
            variant_set['variant_set_id'] = f'variant_set_{variant_set_id:03}'

    def yaml_to_var_list(self):
        try:
            with open(self.par_file) as yaml_file:
                self.config = yaml.full_load(yaml_file)
        except:
            raise Exception("YAML File {} failed to be open".format(self.par_file))
        self.postproc_config_dict()

    def write_to_file(self, sv, bedfile, source_s, source_e, target_s, target_e, transform, event, sv_id):
        assert (not event.symbol.startswith(Symbols.DIS.value))
        if transform == Operations.INS.value:
            transform_length = event.length
        else:
            transform_length = source_e - source_s
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
                       '%d/%d' % (int(sv.hap[0]), int(sv.hap[1])),
                       sv.name,
                       str(sv_id)]
                fout.write("\t".join(row) + "\n")

    @staticmethod
    def symbol_is_inversion(symbol):
        return any(c.islower() for c in symbol)

    @staticmethod
    def export_insertions(chr, start_pos, seq, ins_fasta):
        """
        Exports foreign insertion sequences to separate fasta file, append only
        """
        with open(ins_fasta, "a") as fout_ins:
            fout_ins.write(">{}_{}\n".format(chr, start_pos))
            fout_ins.write("{}\n".format(seq))

    @staticmethod
    def get_event_target_operation(ev, target_events_dict, source_events_dict):
        """
        determines target interval and operation for multi-source events
        """
        # A -> A'
        if ev + Symbols.DUP.value in target_events_dict:
            trg_sym = ev + Symbols.DUP.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), \
                Operations.DUP.value if ev in target_events_dict else Operations.TRA.value
        # A -> a'
        elif ev.lower() + Symbols.DUP.value in target_events_dict:
            trg_sym = ev.lower() + Symbols.DUP.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INVDUP.value
        # A -> a
        elif ev.lower() in target_events_dict:
            trg_sym = ev.lower()
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.INV.value
        # A -> A* (in the case of a custom event in which an event is divergently duplicated)
        elif ev + Symbols.DIV.value in target_events_dict:
            trg_sym = ev + Symbols.DIV.value
            return (target_events_dict[trg_sym].start, target_events_dict[trg_sym].end), Operations.DIV.value
        # A -> A (insertion if source A is undefined, identity otherwise)
        elif ev in target_events_dict:
            return (target_events_dict[ev].start, target_events_dict[ev].end), \
                Operations.INS.value if source_events_dict[ev].start is None else Operations.IDENTITY.value
        # A -> [none]
        elif ev not in [sym[0] for sym in target_events_dict]:
            return (source_events_dict[ev].start, source_events_dict[ev].end), Operations.DEL.value
        # otherwise unknown mapping
        else:
            return (source_events_dict[ev].start, source_events_dict[ev].end), Operations.UNDEFINED.value

    @staticmethod
    def postprocess_record_params(sv, sv_record_info):
        """
        arrange the bed_record parameter dictionaries in order of ascending source interval start position
        and assign order values to the relevant entries
        """
        # for TRA/INS/DUP events with the same target position, 'order' describes the order in which they
        # are compiled (i.e., the order in which they appear in the target sequence)
        order = 0
        ins_pos = None
        for block in sv.target_symbol_blocks:
            for target_event in block:
                if target_event.symbol.startswith(Symbols.DIS.value) or \
                        target_event.symbol in sv_record_info:  # <- prevent collision with A' and A if both in target
                    continue
                src_sym = target_event.symbol[0].upper()
                if sv_record_info[src_sym]['transform'] in NONZERO_ORDER_OPERATIONS:
                    if ins_pos is None:
                        ins_pos = sv_record_info[src_sym]['target_s']
                        order += 1
                    elif sv_record_info[src_sym]['target_s'] == ins_pos:
                        order += 1
                else:
                    ins_pos = None
                    order = 0
                # sv_record_info[src_sym]['order'] = order
        return sorted([params for params in sv_record_info.values()], key=lambda params: params['source_s'])

    @staticmethod
    def get_grammar(sv_config):
        """
        Collect the sorted symbols in corresponding order of length_ranges given
        in variant_sets config entry for a given SV (either predefined or custom)
        """
        if sv_config["type"] != Variant_Type.Custom:
            src_frags, trg_frags = map(utils.reformat_seq, SV_KEY[sv_config["type"]])
            union_frags = set(c[0].upper() if c[0].isalnum() else c for c in src_frags + trg_frags)
        else:
            src_frags, trg_frags = map(utils.reformat_seq, (list(sv_config["source"]), list(sv_config["target"])))
            union_frags = set(c[0].upper() if c[0].isalnum() else c for c in src_frags + trg_frags)
        return sorted([c for c in union_frags if c.isalnum()]) + [c for c in union_frags if not c.isalnum()]

    def export_to_bedpe(self, svs, bedfile, ins_fasta=None, reset_file=True):
        if reset_file:
            utils.reset_file(bedfile)
            if ins_fasta:
                utils.reset_file(ins_fasta)
        for sv_id, sv in enumerate(svs):
            # SVs with multiple source events will be split into multiple bed records (one for each)
            if len(sv.events_dict) == 1:
                ev = list(sv.sv_blocks.target_events_dict.values())[0] if sv.type == Variant_Type.INS\
                        else list(sv.events_dict.values())[0]
                op = self.get_event_target_operation(ev.symbol, sv.sv_blocks.target_events_dict, sv.events_dict)[1]
                record_info = {'source_s': ev.start, 'source_e': ev.end, 'target_s': ev.start, 'target_e': ev.end,
                               'transform': op, 'sv': sv, 'event': ev, 'bedfile': bedfile, 'sv_id': sv_id + 1}
                self.write_to_file(**record_info)
                if op == Operations.INS.value:
                    self.export_insertions(sv.start_chr, ev.start, ev.source_frag, ins_fasta)
            else:
                # multiple source events: source intervals taken from the source events
                # and target intervals taken from corresponding target events (if no match, then deletion)
                sv_record_info = {}
                for ev in sv.events_dict.values():
                    if ev.symbol.startswith(Symbols.DIS.value):
                        continue
                    sv_record_info[ev.symbol] = {'source_s': ev.start, 'source_e': ev.end, 'sv': sv, 'event': ev, 'bedfile': bedfile, 'sv_id': sv_id + 1}
                    (target_s, target_e), operation = self.get_event_target_operation(ev.symbol, sv.sv_blocks.target_events_dict, sv.events_dict)
                    if operation == Operations.INS.value:
                        target_ev = sv.sv_blocks.target_events_dict[ev.symbol]
                        sv_record_info[ev.symbol]['source_s'] = target_s
                        sv_record_info[ev.symbol]['source_e'] = target_e
                        sv_record_info[ev.symbol]['event'] = target_ev
                        self.export_insertions(sv.start_chr, target_s, target_ev.source_frag, ins_fasta)
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
        vcf_file.header.info.add('OVERLAP_EV', number=1, type='String', description="Bool. indicator for the event being"
                                                                                    "placed at an overlap_regions interval")
        vcf_file.header.formats.add('GT', number=1, type='String', description="Genotype")

        vcf_out_file = pysam.VariantFile(vcffile, 'w', header=vcf_file.header)

        for sv in svs:
            zyg = (int(sv.hap[0]), int(sv.hap[1]))
            dispersion_target = None
            if sv.type in DISPERSION_TYPES:
                source_event = sv.events_dict[Symbols.REQUIRED_SOURCE.value]
                disp_event = sv.events_dict['_1']
                rec_start = source_event.start
                rec_end = source_event.end
                if disp_event.start == source_event.end:
                    dispersion_target = disp_event.end
                else:
                    dispersion_target = disp_event.start
            else:
                rec_start = min([frag[1] for frag in sv.changed_fragments])
                rec_end = max(frag[2] for frag in sv.changed_fragments)
            if dispersion_target is not None:
                info_field = {'SVTYPE': sv.type.value, 'SVLEN': rec_end - rec_start, 'TARGET': dispersion_target}
            else:
                if sv.type == Variant_Type.INS:
                    # special case of simple INS: sv length \neq (sv end - sv start)
                    # **pysam will delete END fields that are equal to POS, therefore INS records won't have an END
                    rec_end += 1
                    info_field = {'SVTYPE': sv.type.value, 'SVLEN': sv.events_dict[Symbols.REQUIRED_SOURCE.value].length}
                else:
                    info_field = {'SVTYPE': sv.type.value, 'SVLEN': rec_end - rec_start}
            if sv.overlap_events is not None:
                ovlp_str = '_'.join(ovlp[3] for ovlp in sv.overlap_events.values() if ovlp is not None)
                if ovlp_str != '':
                    info_field['OVERLAP_EV'] = ovlp_str

            vcf_record = vcf_out_file.header.new_record(contig=sv.start_chr, start=rec_start, stop=rec_end,
                                                        alleles=['N', '<%s>' % sv.type.value], id=sv.type.value,
                                                        info=info_field,
                                                        qual=100, filter='PASS',
                                                        samples=[{'GT': zyg}])
            vcf_out_file.write(vcf_record)

        vcf_out_file.close()

    def export_variants_to_fasta(self, id, edits, fasta_out, fasta_file, verbose=False):
        """
        Exports list of changes from simulator to fasta file

        id: chr_id to apply edits to
        edits: list with elements of the form (start, end, new_frag)
        fasta_out: Fasta file to export changes to
        fasta_file: FastaFile with access to reference
        """
        with open(fasta_out, "a") as fout_export:
            if id not in fasta_file.references:
                raise KeyError("ID {} not found in inputted fasta file".format(id))
            if verbose:
                print("New ID: ", id)
            fout_export.write(">" + str(id) + "\n")
            chr_variants = list(edits)
            chr_variants.sort()
            chr_variants.append([fasta_file.get_reference_length(id), fasta_file.get_reference_length(id), ""])
            pos = 0
            for variant in chr_variants:
                var_start, var_end = variant[0], variant[1]
                while pos < var_start:
                    appropriate_buffer = MAX_BUFFER_SIZE if var_start - pos > MAX_BUFFER_SIZE else var_start - pos
                    c = fasta_file.fetch(id, pos, pos + appropriate_buffer)
                    fout_export.write(c)
                    pos += appropriate_buffer
                assert (pos == var_start), "Replacement fragment about to be inserted at position {} instead of var_start {}".format(pos, var_start)
                fout_export.write(variant[2])
                pos = var_end
            fout_export.write("\n")

    def close(self):
        self.fin_export1.close()
        self.fin_export2.close()


def collect_args():
    parser = argparse.ArgumentParser(description='insilicoSV is a software to design and simulate complex structural variants, both novel and known.')
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("-r", "--root", action="store", metavar="DIR", dest="root_dir", help="root directory for all files given")
    parser.add_argument("--random_seed", type=int,
                        help="if non-zero, random seed for random number generation")

    args = parser.parse_args()
    output_dir = os.path.dirname(args.config)
    args_dict = {"config": args.config, "ins_fasta": os.path.join(output_dir, "sim.insertions.fa"),
                 "hap1": os.path.join(output_dir, "sim.hapA.fa"), "hap2": os.path.join(output_dir, "sim.hapB.fa"),
                 "bedpe": os.path.join(output_dir, "sim.bed"), "stats": os.path.join(output_dir, "sim.stats.txt"),
                 "log_file": os.path.join(output_dir, "sim.log"), "random_seed": args.random_seed}
    if args.root_dir:
        for key, curr_path in args_dict.items():
            args_dict[key] = os.path.join(args.root_dir, curr_path)
    return args_dict


