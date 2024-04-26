from collections import defaultdict
from os import write
import logging
import os
import random
import time

from pysam import FastaFile
from pysam import VariantFile

from insilicosv import utils
from insilicosv.processing import FormatterIO, collect_args
from insilicosv.constants import Variant_Type, Zygosity, Symbols
from insilicosv.structural_variant import Structural_Variant

time_start = time.time()


class StatsCollection:
    """collection of information for stats file, if requested"""
    # def __init__(self, chr_lens):
    def __init__(self, chr_ids, chr_lens):
        self.num_heterozygous = 0
        self.num_homozygous = 0
        self.total_svs = 0
        self.active_svs = 0
        self.active_events_chr = dict()
        self.chr_lengths = chr_lens
        self.chr_ids = chr_ids
        self.avg_len = [0, 0]  # Average length of SV events/components
        self.len_frags_chr = dict()  # Lengths of altered fragments within chromosome
        self.sv_types = dict()

        for id in self.chr_ids:
            self.len_frags_chr[id] = 0
            self.active_events_chr[id] = 0

    def get_info(self, svs):
        """
        collects all information for stats file after all edits are completed

        svs: list of Structural_Variant objects
        -> return None
        """
        self.total_svs = len(svs)
        self.min_event_len = min([event.length for sv in svs if sv.active for key, event in sv.events_dict.items() if
                                  not event.symbol.startswith(Symbols.DIS.value)])
        self.max_event_len = max([event.length for sv in svs if sv.active for key, event in sv.events_dict.items() if
                                  not event.symbol.startswith(Symbols.DIS.value)])

        for sv in svs:
            if sv.active:
                self.active_svs += 1
                if sv.name in self.sv_types:
                    self.sv_types[sv.name] += 1
                else:
                    self.sv_types[sv.name] = 1

                if sv.hap[0] and sv.hap[1]:  # homozygous SV
                    self.num_homozygous += 1
                else:  # heterozygous SV
                    self.num_heterozygous += 1

                # add up the lengths of impacted regions on the reference
                for frag in sv.changed_fragments:
                    self.len_frags_chr[frag[0]] += frag[2] - frag[1]

                # count up average length of non-dispersion events
                for symbol in sv.events_dict:
                    if not symbol.startswith(Symbols.DIS.value):
                        event = sv.events_dict[symbol]
                        self.avg_len[0] += event.length
                        self.avg_len[1] += 1
        self.avg_len = self.avg_len[0] // self.avg_len[1] if self.avg_len[1] != 0 else 0

    def export_data(self, fileout):
        """
        Exports all collected data to entered file
        fileout: Location to export stats file
        """
        def write_item(fout, name, item, prefix=""):
            fout.write("{}{}: {}\n".format(prefix, str(name), str(item)))

        os.makedirs(os.path.dirname(fileout), exist_ok=True)
        with open(fileout, "w") as fout:
            fout.write("===== Overview =====\n")
            write_item(fout, "SVs successfully simulated", str(self.active_svs) + "/" + str(self.total_svs))
            for sv_type in self.sv_types:
                write_item(fout, sv_type, self.sv_types[sv_type], prefix="\t- ")
            write_item(fout, "Homozygous SVs", self.num_homozygous)
            write_item(fout, "Heterozygous SVs", self.num_heterozygous)
            write_item(fout, "Average length of SV symbols/components (excluding dispersions)", self.avg_len)
            write_item(fout, "Min length of non-dispersion event", self.min_event_len)
            write_item(fout, "Max length of non-dispersion event", self.max_event_len)
            for id in self.chr_ids:
                fout.write("\n===== {} =====\n".format(id))
                write_item(fout, "Length of sequence", self.chr_lengths[id])
                write_item(fout, "Total impacted length of reference chromosome", self.len_frags_chr[id])


class SV_Simulator:
    def __init__(self, par_file, log_file=None):
        """
        par_file: file location to configuration file (.yaml)
        log_file: location to store log file with diagnostic information if config parameters indicate so
        """
        global time_start
        print("Setting up Simulator...")

        self.formatter = FormatterIO(par_file)
        self.formatter.yaml_to_var_list()
        config = self.formatter.config
        self.ref_file = config['sim_settings']['reference']
        self.ref_fasta = FastaFile(self.ref_file)
        self.svs_config = config['variant_sets']
        self.min_intersv_dist = int(config['sim_settings']['min_intersv_dist']) if \
            'min_intersv_dist' in config['sim_settings'] else 1

        self.sim_settings = config['sim_settings']
        if log_file and "generate_log_file" in self.sim_settings:
            logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG,
                                format='[%(name)s: %(levelname)s - %(asctime)s] %(message)s')
            self.log_to_file("YAML Configuration: {}".format(config))

        # get all chromosome ids
        self.order_ids = self.ref_fasta.references
        self.len_dict = dict()  # stores mapping with key = chromosome, value = chromosome length
        self.filtered_chroms = set()
        for id in self.order_ids:
            chrom_len = self.ref_fasta.get_reference_length(id)
            self.len_dict[id] = chrom_len
            print("Length of chromosome {}: {}".format(id, self.len_dict[id]))
            if 'filter_small_chr' in self.sim_settings and chrom_len < self.sim_settings['filter_small_chr']:
                self.filtered_chroms.add(id)
                print("Filtering chromosome {}: Length of {} below threshold of {}".format(id, chrom_len, self.sim_settings['filter_small_chr']))

        # initialize stats file to be generated after all edits and exporting are finished
        self.stats = StatsCollection(self.order_ids, self.len_dict)

        self.mode = "randomized"
        self.vcf_path = None
        if "vcf_path" in self.svs_config[0]:
            self.mode = "fixed"
            self.vcf_path = self.svs_config[0]["vcf_path"]

        self.svs = []
        self.event_ranges = defaultdict(list)

        self.blacklist_regions = None if "blacklist_regions" not in config.keys() \
            else utils.BlacklistRegions(config)

        self.overlap_regions = None if "overlap_regions" not in config.keys() \
            else utils.OverlapRegions(config, allow_chroms=self.order_ids)

        self.active_svs_total = 0
        self.inactive_svs_total = 0

        self.initialize_svs()

        print("Finished Setting up Simulator in {} seconds\n".format(time.time() - time_start))
        time_start = time.time()

    def __repr__(self):
        return "All structural variants entered into simulator: {}".format(self.svs)

    def log_to_file(self, info, key="DEBUG"):
        # only logs to file if config setting indicates so
        key_to_func = {"DEBUG": logging.debug, "WARNING": logging.warning}
        if "generate_log_file" in self.sim_settings and self.sim_settings["generate_log_file"]:
            key_to_func[key](info)

    def get_rand_chr(self, check_size=None, fixed_chrom=None):
        # random assignment of SV to a chromosome (unless we have a predetermined chromosome for this event)
        valid_chrs = [chrom for chrom, chr_size in self.len_dict.items() if chrom not in self.filtered_chroms and
                      (check_size is None or chr_size >= check_size)]
        if len(valid_chrs) == 0:
            raise Exception("No valid chromosomes remaining")
        rand_id = random.choice(valid_chrs) if fixed_chrom is None else fixed_chrom
        chr_len = self.len_dict[rand_id]

        return rand_id, chr_len

    def extract_event_intervals(self, file_path):
        if file_path.endswith('.vcf') or file_path.endswith('.vcf.gz'):
            for rec in VariantFile(file_path).fetch():
                self.event_ranges[rec.chrom].append((rec.start, rec.stop))
        elif file_path.endswith('.bed'):
            with open(file_path) as f:
                for line in f:
                    chrom, start, stop = tuple(line.split()[:3])
                    self.event_ranges[chrom].append((int(start), int(stop)))
        else:
            raise Exception('extract_event_intervals() called on file of invalid type')

    def process_vcf(self, vcf_path):
        # process vcf containing SVs to be added (deterministically) to reference
        self.active_svs_total = 0
        time_start_local = 0
        vcf = VariantFile(vcf_path)
        for rec in vcf.fetch():
            svtype = Variant_Type(rec.info['SVTYPE']) if 'SVTYPE' in rec.info else Variant_Type(rec.id)
            self.event_ranges[rec.chrom].append((rec.start, rec.stop))
            sv = Structural_Variant(sv_type=svtype, sv_chrom=rec.chrom, mode='fixed', vcf_rec=rec, ref_fasta=self.ref_fasta)
            self.svs.append(sv)
            self.active_svs_total += 1
            self.log_to_file("Intervals {} added to Chromosome \"{}\"".format(self.event_ranges[rec.chrom], rec.chrom))
            time_dif = time.time() - time_start_local
            print("{} SVs successfully placed ========== {} seconds".format(self.active_svs_total, time_dif), end="\r")
            time_start_local = time.time()

    def initialize_svs(self):
        """
        Creates Structural_Variant objects for every SV to simulate and decides zygosity
        self.mode: flag indicating whether SVs are to be randomly generated or read in from VCF
        self.vcf_path: optional path that will be used if mode=="fixed"
        """
        if self.mode == "randomized":
            for sv_config in self.svs_config:
                # collect sv_config-level info
                sv_blacklist_regions = self.blacklist_regions.get_roi(sv_config) \
                    if 'blacklist_region_type' in sv_config else None
                for num in range(sv_config["number"]):
                    repeat_elts = None
                    if 'overlap_type' in sv_config or 'overlap_component' in sv_config:
                        sv_config_identifier = utils.get_sv_config_identifier(sv_config)
                        if sv_config_identifier in self.overlap_regions.svtype_overlap_counts:
                            repeat_elts = self.overlap_regions.get_overlap_intervals(
                                sv_config_identifier, sv_config, partial_overlap=False)
                        elif sv_config_identifier in self.overlap_regions.svtype_partial_overlap_counts:
                            repeat_elts = self.overlap_regions.get_overlap_intervals(
                                sv_config_identifier, sv_config, partial_overlap=True)
                        elif sv_config_identifier in self.overlap_regions.svtype_alu_mediated_counts:
                            repeat_elts = self.overlap_regions.get_alu_mediated_interval(sv_config_identifier)
                        # if SV intended to overlap ROI but no ROI is successfully found, fail to place
                        if repeat_elts is None or set(repeat_elts.values()) == {None}:
                            self.inactive_svs_total += 1
                            continue
                    sv_chrom, chr_len = self.get_rand_chr(fixed_chrom=(
                        None if repeat_elts is None or list(repeat_elts.values()) == [None] else
                        list(repeat_elts.values())[0][0]))
                    if sv_config['type'] == Variant_Type.SNP:
                        sv = Structural_Variant(sv_type=sv_config["type"], sv_chrom=sv_chrom, mode=self.mode, length_ranges=[(1, 1)],
                                                chrom_length=chr_len, blacklist_regions=sv_blacklist_regions)
                    else:
                        sv = Structural_Variant(sv_type=sv_config["type"], sv_chrom=sv_chrom, mode=self.mode,
                                                length_ranges=sv_config["length_ranges"], chrom_length=chr_len,
                                                source=sv_config["source"], target=sv_config["target"],
                                                overlap_events=repeat_elts,
                                                div_prob=(None if 'divergence_prob' not in sv_config else sv_config['divergence_prob']),
                                                blacklist_regions=sv_blacklist_regions)

                    if self.sim_settings.get("homozygous_only", False) or random.randint(0, 1):
                        sv.ishomozygous = Zygosity.HOMOZYGOUS
                        sv.hap = [True, True]
                    else:
                        sv.ishomozygous = Zygosity.HETEROZYGOUS
                        sv.hap = random.choice([[True, False], [False, True]])

                    self.svs.append(sv)
            if not ("prioritize_top" in self.sim_settings and self.sim_settings["prioritize_top"]):
                random.shuffle(self.svs)
        else:  # mode == "fixed"
            self.process_vcf(self.vcf_path)

    def produce_variant_genome(self, fasta1_out, fasta2_out, ins_fasta, bedfile, stats_file=None, initial_reset=True,
                               verbose=False, export_to_file=True):
        """
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        stats_file: whether a stats file summarizing SVs simulated will be generated in same directory the reference genome is located in
        """
        global time_start
        if initial_reset:
            utils.reset_file(fasta1_out)
            utils.reset_file(fasta2_out)
        ref_fasta = self.ref_fasta
        self.apply_transformations(ref_fasta)
        print("Finished SV placements and transformations in {} seconds".format(time.time() - time_start))
        time_start = time.time()
        active_svs = [sv for sv in self.svs if sv.active]
        print("Starting Export Process...")
        for x in range(2):
            edits_dict = dict()
            for id in self.order_ids:
                edits_dict[id] = []
            if x == 0:
                fasta_out = fasta1_out
            elif x == 1:
                fasta_out = fasta2_out
            for sv in active_svs:
                if sv.hap[x]:
                    for frag in sv.changed_fragments:
                        edits_dict[frag[0]].append(frag[1:])
            for id in edits_dict:
                edits_dict[id].sort()
                self.event_ranges[id].sort()
            self.log_to_file("Event Ranges: {}".format(self.event_ranges))
            self.log_to_file("Intervals for hap {}: {}".format(x, edits_dict))
            for id in self.order_ids:
                edits_x = edits_dict[id]
                utils.fail_if_any_overlapping(edits_x)
                self.formatter.export_variants_to_fasta(id, edits_x, fasta_out, ref_fasta, verbose=verbose)
                print("ID {} exported to fasta file {} in {} seconds".format(id, fasta_out, time.time() - time_start))
                time_start = time.time()
        if export_to_file:
            self.formatter.export_to_bedpe(active_svs, bedfile, ins_fasta=ins_fasta, reset_file=initial_reset)
            self.formatter.export_to_vcf(active_svs, self.stats, vcffile=bedfile[:-4]+'.vcf')
        if stats_file:
            self.stats.get_info(self.svs)
            self.stats.export_data(stats_file)

    def choose_rand_pos(self, svs, ref_fasta):
        """
        randomly positions SVs and stores reference fragments in SV events

        svs: list of Structural Variant objects
        ref_fasta: FastaFile with access to reference file
        """
        time_start_local = time.time()
        if self.sim_settings["fail_if_placement_issues"] and self.inactive_svs_total > 0:
            raise Exception("Failure to place SVs at specified ROI-constrained positions")
        for sv in svs:
            tries = 0
            valid = False
            while not valid:
                tries += 1
                valid = True
                if tries > self.sim_settings["max_tries"]:
                    if self.sim_settings["fail_if_placement_issues"]:
                        raise Exception(
                            "Failed to simulate {}, {} / {} SVs successfully simulated (set fail_if_placement_issues "
                            "to False to override placement failures)".format(
                                sv, self.active_svs_total, self.active_svs_total + self.inactive_svs_total))
                    valid = False
                    break
                new_intervals = []
                sv_event_new_vals = []

                # if no overlap, randomly draw start position, otherwise use the overlap to anchor SV placement
                if all(frag.start is None for frag in sv.source_events):
                    start_pos = random.randint(0, sv.chrom_length - sv.req_space)
                else:
                    start_pos = 0
                    for frag in sv.source_events[::-1]:
                        if frag.start is not None:
                            start_pos = frag.start
                        else:
                            start_pos -= frag.length
                sv_start = start_pos
                sv_end = sv_start + sv.req_space
                block_start = sv_start

                for sv_event in sv.source_events:
                    sv_event_source_chr = sv.start_chr
                    sv_event_start, sv_event_end = start_pos, start_pos + sv_event.length
                    frag = ref_fasta.fetch(sv.start_chr, sv_event_start, sv_event_end)
                    sv_event_source_frag = frag
                    sv_event_new_vals.append((sv_event_source_chr, sv_event_start, sv_event_end, sv_event_source_frag))
                    start_pos += sv_event.length

                    if sv_event.symbol.startswith(Symbols.DIS.value):
                        if utils.is_overlapping(self.event_ranges[sv.start_chr],
                                                (block_start - self.min_intersv_dist, sv_event_start + self.min_intersv_dist)) or \
                                (sv.blacklist_region_trees is not None and
                                 sv.blacklist_region_trees[sv.start_chr].overlap(block_start, sv_event_start)):
                            valid = False
                            break
                        new_intervals.append((block_start, sv_event_start))
                        block_start = sv_event_end
                    elif utils.percent_N(frag) > 0.05:
                        valid = False
                        break
                # catches the last (and perhaps only) block in sequence
                if utils.is_overlapping(self.event_ranges[sv.start_chr], (block_start - self.min_intersv_dist,
                                                                          sv_end + self.min_intersv_dist)) or \
                        (sv.blacklist_region_trees is not None and
                         sv.blacklist_region_trees[sv.start_chr].overlap(block_start, sv_end)):
                    valid = False
                    continue
                else:
                    new_intervals.append((block_start, sv_end))

            # adds new SV to simulate only if chosen positions were valid
            if valid:
                self.active_svs_total += 1
                sv.active = True
                sv.start, sv.end = sv_start, sv_end
                for sv_event, sv_event_new_val in zip(sv.source_events, sv_event_new_vals):
                    sv_event.source_chr, sv_event.start, sv_event.end, sv_event.source_frag = sv_event_new_val
                self.log_to_file("Intervals {} added to Chromosome \"{}\" for SV {}".format(new_intervals, sv.start_chr, sv))
                self.event_ranges[sv.start_chr].extend(new_intervals)
                sv.assign_locations(sv.start)
            else:
                self.inactive_svs_total += 1
                if tries != self.sim_settings["max_tries"] + 1:
                    self.log_to_file("{} only got {} tries instead of the max {}".format(sv, tries, self.sim_settings[
                        "max_tries"] + 1), key="WARNING")

            time_dif = time.time() - time_start_local
            print(
                "{} / {} SVs successfully placed ========== {} / {} SVs unsuccessfully placed, {} tries, {} seconds".format(
                    self.active_svs_total, self.active_svs_total + self.inactive_svs_total,
                    self.inactive_svs_total, self.active_svs_total + self.inactive_svs_total, tries, time_dif), end="\r")
            time_start_local = time.time()

    def apply_transformations(self, ref_fasta):
        """
        Randomly chooses positions for all SVs and carries out all edits
        Populates event classes within SVs with reference fragments and start & end positions
        Stores list of changes, which each have an interval and a sequence to substitute the reference frag with, in SV

        ref_fasta: FastaFile with access to reference
        mode: flag indicating whether we're adding SVs to the reference in a randomized or deterministic way
        """
        if self.mode == "randomized":
            # select random positions for SVs
            self.choose_rand_pos(self.svs, ref_fasta)
            print()

        total = 0
        for sv in self.svs:
            if sv.active:
                sv.change_fragment()
                total += 1
                self.log_to_file("Events Dict after all edits: {} ".format(sv.events_dict))

    def close(self):
        self.ref_fasta.close()

def run_insilicosv():
    sim_start = time.time()

    args = collect_args()
    if args["random_seed"]:
        random.seed(args["random_seed"])
        print(f"Random seed: {args['random_seed']}")
    yaml_in = args["config"]
    fasta1_out = args["hap1"]
    fasta2_out = args["hap2"]
    ins_fasta = args["ins_fasta"]
    bed_out = args["bedpe"]
    stats_file = args["stats"]
    log_file = args["log_file"]

    sim = SV_Simulator(args["config"], log_file=args["log_file"])
    sim.produce_variant_genome(args["hap1"], args["hap2"], args["ins_fasta"], args["bedpe"], args["stats"],
                               verbose=False)
    print("Simulation completed in {} seconds".format(time.time() - sim_start))

if __name__ == "__main__":
    run_insilicosv()
