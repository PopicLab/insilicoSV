from os import write
import random
from pysam import FastaFile
from processing import FormatterIO, collect_args
import utils
from constants import *
from structural_variant import Structural_Variant, Event
import tracemalloc  # only for testing
import sys
import logging
import time

time_start = time.time()

class StatsCollection():
    '''collection of information for stats file, if requested'''
    def __init__(self, chr_ids, chr_lens):
        self.num_heterozygous = 0
        self.num_homozygous = 0
        self.total_svs = 0
        self.active_svs = 0
        self.active_events_chr = dict()
        self.chr_ids = chr_ids
        self.chr_lengths = chr_lens
        self.avg_len = [0,0]       # Average length of SV events/components
        self.len_frags_chr = dict()  # Lengths of altered fragments within chromosome
        self.sv_types = dict()

        for id in self.chr_ids:
            self.len_frags_chr[id] = 0
            self.active_events_chr[id] = 0
    
    def get_info(self, svs):
        '''
        collects all information for stats file after all edits are completed

        svs: list of Structural_Variant objects
        -> return None
        '''

        self.total_svs = len(svs)
        self.min_event_len = min([event.length for sv in svs if sv.active for key, event in sv.events_dict.items() if not event.symbol.startswith(Symbols.DIS.value)])
        self.max_event_len = max([event.length for sv in svs if sv.active for key, event in sv.events_dict.items() if not event.symbol.startswith(Symbols.DIS.value)])

        for sv in svs:
            if sv.active: # only collect information for SVs that were successfully simulated
                self.active_svs += 1
                # count up SVs of each type
                if sv.name in self.sv_types:
                    self.sv_types[sv.name] += 1
                else:
                    self.sv_types[sv.name] = 1

                # count up zygosity
                if sv.hap[0] and sv.hap[1]: # homozygous SV
                    self.num_homozygous += 1
                else:  # heterozygous SV
                    self.num_heterozygous += 1

                # add up the lengths of impacted regions on the reference
                for frag in sv.changed_fragments:
                    self.len_frags_chr[frag[0]] += frag[2] - frag[1]    # frag[0] = chromosome id, frag[1] = start position, frag[2] = end position of edit
                    
                # count up average length of non-dispersion events
                for symbol in sv.events_dict:
                    if not symbol.startswith(Symbols.DIS.value):   
                        event = sv.events_dict[symbol]
                        self.avg_len[0] += event.length
                        self.avg_len[1] += 1
        if self.avg_len[1] != 0:   # avoid ZeroDivisionError
            self.avg_len = self.avg_len[0] // self.avg_len[1]
        else:
            self.avg_len = 0
    def export_data(self, fileout):
        '''
        Exports all collected data to entered file
        fileout: Location to export stats file
        '''
        def write_item(fout, name, item, prefix = ""):
            fout.write("{}{}: {}\n".format(prefix, str(name), str(item)))
        with open(fileout, "w") as fout:
            fout.write("===== Overview =====\n")
            write_item(fout, "SVs successfully simulated", str(self.active_svs) + "/" + str(self.total_svs))
            for sv_type in self.sv_types:
                write_item(fout, sv_type, self.sv_types[sv_type], prefix = "\t- ")
            write_item(fout, "Homozygous SVs", self.num_homozygous)
            write_item(fout, "Heterozygous SVs", self.num_heterozygous)
            write_item(fout, "Average length of SV symbols/components (excluding dispersions)", self.avg_len)
            write_item(fout, "Min length of non-dispersion event", self.min_event_len)
            write_item(fout, "Max length of non-dispersion event", self.max_event_len)
            for id in self.chr_ids:
                fout.write("\n===== {} =====\n".format(id))
                write_item(fout, "Length of sequence", self.chr_lengths[id])
                write_item(fout, "Total impacted length of reference chromosome", self.len_frags_chr[id])

class SV_Simulator():
    def __init__(self, ref_file, par_file, random_gen = random, log_file=None):
        '''
        ref_file: file location to reference genome (.fasta or .fna or .fa)
        par_file: file location to configuration file (.yaml)
        random_gen: only relevant for unittesting
        log_file: location to store log file with diagnostic information if config parameters indicate so
        '''
        global time_start
        print("Setting up Simulator...")

        # create FastaFile with access to the reference file
        self.ref_file = ref_file
        #utils.remove_file(self.ref_file + ".fai")    # if ref file was changed, then previously generated index file will not be updated
        self.ref_fasta = FastaFile(ref_file)        # FastaFile generates a new index file if one is not found

        # get all chromosome ids
        self.order_ids = self.ref_fasta.references
        self.len_dict = dict()  # stores mapping with key = chromosome, value = chromosome length
        for id in self.order_ids:
            self.len_dict[id] = self.ref_fasta.get_reference_length(id)
            print("Length of chromosome {}: {}".format(id, self.len_dict[id]))
        
        # initialize stats file to be generated after all edits and exporting are finished
        self.stats = StatsCollection(self.order_ids, self.len_dict)

        # process config file for SVs
        self.formatter = FormatterIO(par_file)
        config = self.formatter.yaml_to_var_list()
        self.svs_config = config.SVs    # config for what SVs to simulate
        self.sim_settings = config.sim_settings
        if log_file and self.sim_settings["generate_log_file"]:
            logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG, format='[%(name)s: %(levelname)s - %(asctime)s] %(message)s')
            self.log_to_file("YAML Configuration: {}".format(config.__dict__))

        # create all SVs
        self.svs = []
        self.initialize_svs(random_gen = random_gen)

        print("Finished Setting up Simulator in {} seconds\n".format(time.time() - time_start))
        time_start = time.time()
    
    def __repr__(self):
        return "All structural variants entered into simulator: {}".format(self.svs)
    
    def log_to_file(self, info, key = "DEBUG"):
        # only logs to file if config setting indicates so
        log_func = [logging.debug, logging.warning]
        key_to_func = {"DEBUG": logging.debug, "WARNING": logging.warning}
        if self.sim_settings["generate_log_file"]:
            key_to_func[key](info)
            #logging.debug(info)

    def initialize_svs(self, random_gen = random):
        '''
        Creates Structural_Variant objects for every SV to simulate and decides zygosity
        '''
        for sv_config in self.svs_config:
            for num in range(sv_config["number"]):
                sv = Structural_Variant(sv_config["type"], sv_config["length_ranges"], source=sv_config["source"], target=sv_config["target"]) # inputs: SV type, range of lengths
                draw = random_gen.randint(1,3)
                
                if draw == 3:   # sv applies to both haplotypes
                    sv.ishomozygous = Zygosity.HOMOZYGOUS
                    sv.hap = [True, True]
                elif draw == 2:
                    sv.ishomozygous = Zygosity.HETEROZYGOUS
                    sv.hap = [False, True]
                elif draw == 1:
                    sv.ishomozygous = Zygosity.HETEROZYGOUS
                    sv.hap = [True, False]
                self.svs.append(sv) 
        # shuffle svs if we are not prioritizing the simulation of the SVs inputted first
        if not self.sim_settings["prioritize_top"]:
            random.shuffle(self.svs)

    def produce_variant_genome(self, fasta1_out, fasta2_out, ins_fasta, bedfile, stats_file = None, initial_reset = True, random_gen = random, verbose = False):
        '''
        initial_reset: boolean to indicate if output file should be overwritten (True) or appended to (False)
        random_gen: only relevant in testing
        stats_file: whether a stats file summarizing SVs simulated will be generated in same directory the reference genome is located in
        '''
        global time_start
        if initial_reset:
            utils.reset_file(fasta1_out)
            utils.reset_file(fasta2_out)

        # edit chromosome
        ref_fasta = self.ref_fasta
        self.apply_transformations(ref_fasta, random_gen)
        print("Finished SV placements and transformations in {} seconds".format(time.time() - time_start))
        time_start = time.time()

        # organize edits and export
        active_svs = [sv for sv in self.svs if sv.active]
        
        print("Starting Export Process...")
        for x in range(2):
            edits_dict = dict()
            for id in self.order_ids:
                edits_dict[id] = []

            if x == 0:
                fasta_out = fasta1_out  # where to write new genome to
            elif x == 1:
                fasta_out = fasta2_out
            for sv in active_svs:
                if sv.hap[x]:
                    for frag in sv.changed_fragments:
                        #print("Frag: {}, SV: {}".format(frag[0:3], sv))
                        edits_dict[frag[0]].append(frag[1:])
            
            # only for debugging
            for id in edits_dict:
                edits_dict[id].sort()
                self.event_ranges[id].sort()
            self.log_to_file("Event Ranges: {}".format(self.event_ranges))
            self.log_to_file("Intervals for hap {}: {}".format(x, edits_dict))

            for id in self.order_ids:
                # account for homozygous and heterogeneous variants
                edits_x = edits_dict[id]
                #print("Edits_x: ", edits_x)
                utils.fail_if_any_overlapping(edits_x)
                
                # export edited chromosomes to FASTA files
                self.formatter.export_variants_to_fasta(id, edits_x, fasta_out, ref_fasta, verbose = verbose)

                print("ID {} exported to fasta file {} in {} seconds".format(id, fasta_out, time.time() - time_start))
                time_start = time.time()

        # export variant data to BED file
        self.formatter.export_to_bedpe(active_svs, bedfile, ins_fasta, reset_file = initial_reset)
        initial_reset = False

        # create and export stats file
        if stats_file:
            self.stats.get_info(self.svs)
            self.stats.export_data(stats_file)
        return True
    
    def choose_rand_pos(self, svs, ref_fasta, random_gen = random, verbose = False):
        '''
        randomly positions SVs and stores reference fragments in SV events

        svs: list of Structural Variant objects
        ref_fasta: FastaFile with access to reference file
        random_gen: only relevant for unittesting
        -> returns list of tuples, represents position ranges for non-dispersion events
        '''

        def percent_N(seq):
            total = 0
            if len(seq) == 0:  # avoid ZeroDivisionError
                return 0
            for char in seq:
                if char == "N":
                    total += 1
            return total / len(seq)
        def generate_seq(length):
            # helper function for insertions
            # generates random sequence of bases of given length
            rand_seq = ""
            base_map = {1:"A", 2: "T", 3: "G", 4: "C"}
            for x in range(length):
                rand_seq += base_map[random_gen.randint(1,4)]
            return rand_seq
        def get_rand_chr(check_size=None):
            # random assignment of SV to a chromosome
            # -> returns random chromosome and its length
            valid_chrs = self.order_ids
            if check_size != None and self.sim_settings["filter_small_chr"]:   # if parameter setting on, only choose random chr from those that are big enough
                valid_chrs = [chr for chr, chr_size in self.len_dict.items() if chr_size > check_size]
            if len(valid_chrs) == 0:
                raise Exception("SVs are too big for the reference!")

            rand_id = valid_chrs[random_gen.randint(0, len(valid_chrs)-1)]
            chr_len = self.len_dict[rand_id]
            chr_event_ranges = self.event_ranges[rand_id]
            assert rand_id != None

            return rand_id, chr_len, chr_event_ranges
        
        # maintain separate event ranges for different chromosomes
        self.event_ranges = dict()
        for id in self.order_ids:
            self.event_ranges[id] = []

        active_svs_total = 0
        inactive_svs_total = 0
        time_start_local = 0
        for sv in svs:
            tries = 0 # number of attempts to place sv randomly
            valid = False
            while not valid:
                tries += 1
                valid = True
                
                if tries > self.sim_settings["max_tries"]:
                    if self.sim_settings["fail_if_placement_issues"]:  # user can set a setting to fail if a single SV was not able to be positioned
                        raise Exception("Failed to simulate {}, {} / {} SVs successfully simulated (set fail_if_placement_issues to False to override placement failures)".format(sv, active_svs_total, len(svs)))
                    #print("Failure to simulate \"{}\"".format(sv))
                    valid = False
                    break
                
                rand_id, chr_len, chr_event_ranges = get_rand_chr(check_size=sv.req_space)
                if chr_len - sv.req_space - 1 <= 0:  # SV is too big for current chromosome
                    valid = False
                    continue
                else:
                    start_pos = random_gen.randint(0, chr_len - sv.req_space - 1)
                    # define the space in which SV operates
                    # we now also know where the target positions lie since we know the order and length of events
                    new_intervals = []  # tracks new ranges of blocks
                    sv.start, sv.start_chr = start_pos, rand_id
                    sv.end = sv.start + sv.req_space
                    block_start = sv.start

                    for sv_event in sv.source_events:

                        # store start and end position and reference fragment
                        sv_event.start, sv_event.end = start_pos, start_pos + sv_event.length
                        sv_event.source_chr = rand_id
                        frag = ref_fasta.fetch(rand_id, sv_event.start, sv_event.end)
                        sv_event.source_frag = frag
                        start_pos += sv_event.length

                        # dispersion event should not impact whether position is valid or not, given that spacing is already guaranteed
                        if sv_event.symbol.startswith(Symbols.DIS.value):
                            # check to see if chosen spot is a valid position
                            #print("is_overlapping({}, {}) -> {}".format(chr_event_ranges, (block_start, sv_event.start), utils.is_overlapping(chr_event_ranges, (block_start, sv_event.start))))
                            if utils.is_overlapping(chr_event_ranges, (block_start, sv_event.start)):   # sv_event.start is the end of the current block
                                valid = False
                                break     # if ANY of the non-dispersion events within SV are in an invalid position, then immediately fail the try
                            new_intervals.append((block_start, sv_event.start))
                            block_start = sv_event.end   # remember that the end is actually the position right after the sv_event
                        elif percent_N(frag) > 0.05:
                            valid = False
                            break
                    # catches the last (and perhaps only) block in sequence
                    if utils.is_overlapping(chr_event_ranges, (block_start, sv.end)):  
                        valid = False
                        continue
                    else:
                        new_intervals.append((block_start, sv.end))

            # adds new SV to simulate only if chosen positions were valid
            if valid:
                active_svs_total += 1
                sv.active = True
                self.log_to_file("Intervals {} added to Chromosome \"{}\" for SV {}".format(new_intervals, rand_id, sv))
                chr_event_ranges.extend(new_intervals)

                # populates insertions with random sequence - these event symbols only show up in target transformation
                for event in sv.events_dict.values():
                    if event.source_frag == None and event.length > 0:
                        event.source_frag = generate_seq(event.length)
            else:
                inactive_svs_total += 1
                if tries != self.sim_settings["max_tries"]+1:
                    self.log_to_file("{} only got {} tries instead of the max {}".format(sv, tries, self.sim_settings["max_tries"]+1), key = "WARNING")

            time_dif = time.time() - time_start_local
            print("{} / {} SVs successfully placed ========== {} / {} SVs unsuccessfully placed, {} tries, {} seconds".format(active_svs_total, len(svs), inactive_svs_total, len(svs), tries, time_dif), end = "\r")
            time_start_local = time.time()
        
        print("\n")
        return self.event_ranges

    def apply_transformations(self, ref_fasta, random_gen = random):
        '''
        Randomly chooses positions for all SVs and carries out all edits
        Populates event classes within SVs with reference fragments and start & end positions
        Stores list of changes, which each have an interval and a sequence to substitute the reference frag with, in SV
        
        ref_fasta: FastaFile with access to reference
        '''
        # select random positions for SVs
        self.choose_rand_pos(self.svs, ref_fasta, random_gen)

        print("Starting edit process...")
        active_svs_total = sum([1 for sv in self.svs if sv.active])
        total = 0
        for sv in self.svs:
            if sv.active:
                # make edits and store in sv object
                sv.change_fragment()
                total += 1
                print("{} / {} SVs successfully edited".format(total, active_svs_total), end = "\r")
                self.log_to_file("Events Dict after all edits: {} ".format(sv.events_dict))
    
    def close(self):
        self.ref_fasta.close()

if __name__ == "__main__":
    sim_start = time.time()

    #tracemalloc.start()

    args = collect_args()
    fasta_in = args["ref"]
    yaml_in = args["config"]
    fasta1_out = args["hap1"]
    fasta2_out = args["hap2"]
    ins_fasta = args["ins_fasta"]
    bed_out = args["bedpe"]
    stats_file = args["stats"]
    log_file = args["log_file"]

    '''fasta_in = "debugging/inputs/test.fna"
    yaml_in = "debugging/inputs/par.yaml"
    fasta1_out = "debugging/inputs/test1_out.fna"
    fasta2_out = "debugging/inputs/test2_out.fna"
    bed_out = "debugging/inputs/out.bed"
    stats_file = "debugging/inputs/stats.txt"'''

    sim = SV_Simulator(args["ref"], args["config"], log_file=args["log_file"])
    sim.produce_variant_genome(args["hap1"], args["hap2"], args["ins_fasta"], args["bedpe"], args["stats"], verbose = False)
    #print(str(sim))
    print("Simulation completed in {} seconds".format(time.time() - sim_start))

    #current, peak = tracemalloc.get_traced_memory()
    #print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    #tracemalloc.stop()


