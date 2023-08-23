import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from simulate import SV_Simulator
from pysam import FastaFile
import yaml
import utils


class TestObject():
    def __init__(self, ref, par, hap1, hap2, bed):
        # ref: list, [file location, fasta_contents_as_dictionary]
        # par: list, [file location, yaml_contents_as_dictionary]
        # hap1: str, location to first output fasta file
        # hap2: str, location to second output fasta file
        # bed: str, location to BEDPE file
        self.ref = ref[0]
        self.ref_content = ref[1]
        self.par = par[0]
        self.par_content = par[1]
        self.hap1 = hap1
        self.hap2 = hap2
        self.bed = bed

    def initialize_files(self):
        # initialize all files, delete old test files
        self.remove_test_files()
        self.write_fasta_file(self.ref_content)
        self.write_yaml_file(self.par_content)

    def write_fasta_file(self, chr_to_sequence):
        with open(self.ref, "w") as fout:
            for chr in chr_to_sequence:
                fout.write(">" + chr + "\n")
                fout.write(chr_to_sequence[chr] + "\n")

    def write_yaml_file(self, config_as_dict):
        with open(self.par, "w") as fout:
            yaml.dump(config_as_dict, fout, default_flow_style=False)

    def remove_test_files(self):
        # remove all index files and output files so they can be regenerated
        # if index file not removed, then they won't be remade even if fasta file is altered

        files = [self.ref, self.ref + ".fai", self.par, self.hap1, self.hap1 + ".fai", self.hap2, self.hap2 + ".fai",
                 self.bed]  # remove reference's index (.fai) file because new one should be generated
        for file in files:
            utils.remove_file(file)

    def get_actual_frag(self, return_haps='hap1'):
        # return_haps: ['hap1', 'hap2', 'both'] for control over which haplotype we check
        fasta_out_1 = FastaFile(self.hap1)  # also generates a .fai file
        fasta_out_2 = FastaFile(self.hap2)
        if return_haps == 'hap1':
            return fasta_out_1.fetch(fasta_out_1.references[0], 0,
                                     fasta_out_1.get_reference_length(fasta_out_1.references[0]))
        elif return_haps == 'hap2':
            return fasta_out_2.fetch(fasta_out_2.references[0], 0,
                                     fasta_out_2.get_reference_length(fasta_out_2.references[0]))
        else:
            return (fasta_out_1.fetch(fasta_out_1.references[0], 0,
                                      fasta_out_1.get_reference_length(fasta_out_1.references[0])),
                    fasta_out_2.fetch(fasta_out_2.references[0], 0,
                                      fasta_out_2.get_reference_length(fasta_out_2.references[0])))


class TestSVSimulator(unittest.TestCase):
    def setUp(self):

        # runs before every test
        ref_file = "test/inputs/test.fna"
        par = "test/inputs/par.yaml"

        hap1 = "test/inputs/test1.fna"
        hap2 = "test/inputs/test2.fna"
        bed = "test/inputs/out.bed"

        test_overlap_bed = "test/inputs/example_overlap_events.bed"
        test_overlap_bed_2 = "test/inputs/example_overlap_events_2.bed"

        self.test_objects_no_dis = [TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINVdup", "number": 1, "max_length": 5,
                                                    "min_length": 5}]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINVdup", "number": 1, "max_length": 5, "min_length": 5},
                                                   {"type": "delINVdel", "number": 1, "min_length": 5, "max_length": 5},
                                                   {"type": "dupINVdup", "number": 1, "min_length": 5, "max_length": 5}
                                                   ]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGAGTCAGGGAGCAAAAAAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "dupINVdel", "number": 1, "max_length": 5, "min_length": 5},
                                                   {"type": "delINV", "number": 1, "min_length": 5, "max_length": 5},
                                                   {"type": "INVdel", "number": 1, "min_length": 5, "max_length": 5}
                                                   ]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "dup_INV", "number": 1, "max_length": 5, "min_length": 5},
                                                   {"type": "INV_dup", "number": 1, "min_length": 5, "max_length": 5}
                                                   ]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINVdup", "number": 1, "max_length": 5, "min_length": 5},
                                                   {"avoid_intervals": "test/inputs/example_avoid_interval.vcf"}]}],
                                               hap1, hap2, bed),
                                    # small ref for testing three-part events
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "dupINVdup",
                                                    "number": 1,
                                                    "min_length": {2, 2, 2},
                                                    "max_length": {2, 2, 2}}]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINVdel",
                                                    "number": 1,
                                                    "min_length": {2, 2, 2},
                                                    "max_length": {2, 2, 2}}]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINVdup",
                                                    "number": 1,
                                                    "min_length": {2, 2, 2},
                                                    "max_length": {2, 2, 2}}]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "dupINVdel",
                                                    "number": 1,
                                                    "min_length": {2, 2, 2},
                                                    "max_length": {2, 2, 2}}]}],
                                               hap1, hap2, bed),
                                    # objects for delINV and INVdel
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "delINV",
                                                    "number": 1,
                                                    "min_length": {3, 3},
                                                    "max_length": {3, 3}}]}],
                                               hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "INVdel",
                                                    "number": 1,
                                                    "min_length": {3, 3},
                                                    "max_length": {3, 3}}]}],
                                               hap1, hap2, bed),
                                    # object for inverted duplication
                                    TestObject([ref_file, {"Chromosome19": "CGT"}],
                                               [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                   {"type": "INVdup",
                                                    "number": 1,
                                                    "min_length": 3,
                                                    "max_length": 3}]}],
                                               hap1, hap2, bed)
                                    ]
        # test objects for bidirectional tests
        self.test_dispersion_objects = [TestObject([ref_file, {"Chromosome19": "CT"}],
                                                   [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                       {"type": "TRA",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   hap1, hap2, bed),
                                        TestObject([ref_file, {"Chromosome19": "CT"}],
                                                   [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                       {"type": "dDUP",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   hap1, hap2, bed),
                                        TestObject([ref_file, {"Chromosome19": "CT"}],
                                                   [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                       {"type": "INV_dDUP",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   hap1, hap2, bed),
                                        TestObject([ref_file, {"Chromosome19": "CTG"}],
                                                   [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                       {"type": "dDUP_iDEL",
                                                        "number": 1,
                                                        "min_length": [1, 1, 1],
                                                        "max_length": [1, 1, 1]}]}],
                                                   hap1, hap2, bed),
                                        TestObject([ref_file, {"Chromosome19": "CTG"}],
                                                   [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                       {"type": "INS_iDEL",
                                                        "number": 1,
                                                        "min_length": [1, 1, 1],
                                                        "max_length": [1, 1, 1]}]}],
                                                   hap1, hap2, bed)
                                        ]
        self.test_objects_ins = [TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                            [par, {"sim_settings": {"prioritize_top": True}, "SVs": [
                                                {"type": "INS", "number": 1, "min_length": 5, "max_length": 5},
                                                {"type": "delINV", "number": 1, "min_length": 5, "max_length": 5},
                                                {"type": "INS", "number": 1, "min_length": 5, "max_length": 5}]}],
                                            hap1, hap2, bed)]

        # --------- simple event test objects -----------
        self.test_objects_simple_dels = [TestObject([ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "DEL", "number": 1,
                                                                    "min_length": 13, "max_length": 13}]}],
                                                    hap1, hap2, bed),
                                         TestObject([ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "DEL", "number": 1,
                                                                    "min_length": 14, "max_length": 14}]}],
                                                    hap1, hap2, bed)]

        self.test_objects_simple_dups = [TestObject([ref_file, {"Chromosome19": "CA"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "DUP", "number": 1,
                                                                    "min_length": 2, "max_length": 2}]}],
                                                    hap1, hap2, bed),
                                         TestObject([ref_file, {"Chromosome19": "CAT"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "DUP", "number": 1,
                                                                    "min_length": 2, "max_length": 2}]}],
                                                    hap1, hap2, bed),
                                         TestObject([ref_file, {"Chromosome19": "C"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "DUP", "number": 1,
                                                                    "min_length": 1, "max_length": 1}]}],
                                                    hap1, hap2, bed)]

        self.test_objects_simple_inss = [TestObject([ref_file, {"Chromosome19": "CA"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "INS", "number": 1,
                                                                    "min_length": 5, "max_length": 5}]}],
                                                    hap1, hap2, bed)]

        self.test_objects_simple_invs = [TestObject([ref_file, {"Chromosome19": "CA"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "INV", "number": 1,
                                                                    "min_length": 2, "max_length": 2}]}],
                                                    hap1, hap2, bed),
                                         TestObject([ref_file, {"Chromosome19": "C"}],
                                                    [par, {"sim_settings": {"prioritize_top": True},
                                                           "SVs": [{"type": "INV", "number": 1,
                                                                    "min_length": 1, "max_length": 1}]}],
                                                    hap1, hap2, bed)]
        # ---------- test objects for overlap-aware event placement ------------
        self.test_objects_overlap_simple = [TestObject([ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {"bed": test_overlap_bed},
                                                              "SVs": [{"type": "DEL", "number": 1,
                                                                       "min_length": 2, "max_length": 2,
                                                                       "num_overlap": 1}]}],
                                                       hap1, hap2, bed),
                                            TestObject([ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {"bed": test_overlap_bed},
                                                              "SVs": [{"type": "DUP", "number": 1,
                                                                       "min_length": 2, "max_length": 2,
                                                                       "num_overlap": 1},
                                                                      {"type": "INV", "number": 1,
                                                                       "min_length": 2, "max_length": 2,
                                                                       "num_overlap": 1}]}],
                                                       hap1, hap2, bed),
                                            # combine two input files, filter all but one event by type
                                            TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {"bed": [test_overlap_bed, test_overlap_bed_2],
                                                                                 "allow_types": ["L1PA15"]},
                                                              "SVs": [{"type": "DEL", "number": 1,
                                                                       "min_length": 2, "max_length": 5,
                                                                       "num_overlap": 1}]}],
                                                       hap1, hap2, bed),
                                            # combine two input files, filter all by length
                                            TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {"bed": [test_overlap_bed, test_overlap_bed_2]},
                                                              "SVs": [{"type": "DEL", "number": 1,
                                                                       "min_length": 10, "max_length": 10,
                                                                       "num_overlap": 1}]}],
                                                       hap1, hap2, bed),
                                            # combine two input files, filter all by chromosome
                                            TestObject([ref_file, {"chr19": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [test_overlap_bed, test_overlap_bed_2]},
                                                              "SVs": [{"type": "DEL", "number": 1,
                                                                       "min_length": 10, "max_length": 10,
                                                                       "num_overlap": 1}]}],
                                                       hap1, hap2, bed),
                                            # type-specific num_overlap params
                                            TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [test_overlap_bed, test_overlap_bed_2],
                                                                  "allow_types": ["L1HS", "ALR/Alpha"]},
                                                              "SVs": [{"type": "DEL", "number": 4,
                                                                       "min_length": 1, "max_length": 10,
                                                                       "num_overlap": [2, 1]}]}],
                                                       hap1, hap2, bed),
                                            # type-specific num_overlap param > num available (ALR)
                                            TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                       [par, {"sim_settings": {"prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [test_overlap_bed, test_overlap_bed_2],
                                                                  "allow_types": ["L1HS", "ALR/Alpha"]},
                                                              "SVs": [{"type": "DEL", "number": 5,
                                                                       "min_length": 1, "max_length": 5,
                                                                       "num_overlap": [2, 3]}]}],
                                                       hap1, hap2, bed)
                                            ]
        self.test_objects_overlap_cplx = [TestObject([ref_file, {"chr21": "CTGAT"}],
                                                     [par, {"sim_settings": {"prioritize_top": True},
                                                            "overlap_events": {"bed": [test_overlap_bed, test_overlap_bed_2],
                                                                               "allow_types": ["L1HS"]},
                                                            "SVs": [{"type": "dDUP", "number": 1,
                                                                     "min_length": [2, 1],
                                                                     "max_length": [2, 1],
                                                                     "num_overlap": 1}]}],
                                                     hap1, hap2, bed),
                                          TestObject([ref_file, {"chr21": "CTGATATGGAC"}],
                                                     [par, {"sim_settings": {"prioritize_top": True},
                                                            "overlap_events": {
                                                                "bed": [test_overlap_bed, test_overlap_bed_2],
                                                                "allow_types": ["L1HS", "AluSz6"]},
                                                            "SVs": [{"type": "TRA", "number": 1,
                                                                     "min_length": [4, 1],
                                                                     "max_length": [6, 1],
                                                                     "num_overlap": 1},
                                                                    {"type": "INV_dDUP", "number": 1,
                                                                     "min_length": [1, 1],
                                                                     "max_length": [1, 1],
                                                                     "num_overlap": 1}
                                                                    ]}],
                                                     hap1, hap2, bed)
                                          ]

        # ---------- test objects for divergence event ------------
        self.test_objects_divergence_event = [TestObject([ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [par, {"sim_settings": {"prioritize_top": True},
                                                                "SVs": [{"type": "DIVERGENCE", "number": 1,
                                                                         "min_length": 5, "max_length": 5}]}],
                                                         hap1, hap2, bed)]

        self.test_objects_filter_chroms = [TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                                                  "chr20": "CTCCGT"}],
                                                      [par, {"sim_settings": {"prioritize_top": True,
                                                                              "filter_small_chr": 10},
                                                             "SVs": [{"type": "DEL", "number": 1,
                                                                      "max_length": 3,
                                                                      "min_length": 3}]}],
                                                      hap1, hap2, bed),
                                           TestObject([ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                                                  "chr20": "CTCCGT"}],
                                                      [par, {"sim_settings": {"prioritize_top": True,
                                                                              "filter_small_chr": 50},
                                                             "SVs": [{"type": "DEL", "number": 1,
                                                                      "max_length": 3,
                                                                      "min_length": 3}]}],
                                                      hap1, hap2, bed)]

        self.test_objects_req_space = [TestObject([ref_file, {"chr21": "CTCCGT"}],
                                                  [par, {"sim_settings": {"prioritize_top": True},
                                                         "SVs": [{"type": "DEL", "number": 1,
                                                                  "max_length": 9,
                                                                  "min_length": 9}]}],
                                                  hap1, hap2, bed)]

    def test_is_overlapping(self):
        # non-insertion cases
        self.assertFalse(utils.is_overlapping([(3, 4), (5, 10)], (4, 5)))
        self.assertTrue(utils.is_overlapping([(3, 4), (5, 10)], (3, 5)))
        # insertion cases
        self.assertFalse(utils.is_overlapping([(3, 4), (5, 10)], (4, 4)))
        self.assertFalse(utils.is_overlapping([(3, 4), (5, 10), (10, 15)], (10, 10)))
        self.assertFalse(utils.is_overlapping([(3, 4), (5, 10)], (10, 10)))
        self.assertFalse(utils.is_overlapping([(3, 4), (20, 20), (5, 10)], (20, 20)))
        self.assertFalse(utils.is_overlapping([(3, 4), (20, 20), (5, 10)], (20, 21)))
        self.assertFalse(utils.is_overlapping([(3, 4), (20, 20), (5, 10)], (19, 20)))
        self.assertFalse(utils.is_overlapping([(3, 4), (20, 20), (5, 10)], (21, 21)))
        self.assertFalse(utils.is_overlapping([(3, 4), (5, 10)], (5, 5)))

    # helper method for tests where the output will be in a known list of possibilities
    def helper_test_known_output_svs(self, config_event_obj, target_frags=None):
        # target_frags: optional input frags to be checked for match with output frags
        config = config_event_obj
        config.initialize_files()
        curr_sim = SV_Simulator(config.ref, config.par)
        curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
        changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
        config.remove_test_files()
        if target_frags is not None:
            self.assertTrue(changed_frag_1 in target_frags or changed_frag_2 in target_frags)
        return changed_frag_1, changed_frag_2

    def test_simple_deletions(self):
        self.helper_test_known_output_svs(self.test_objects_simple_dels[0], ['CA', 'CT', 'AT'])
        self.helper_test_known_output_svs(self.test_objects_simple_dels[1], ['C', 'T'])

    def test_simple_duplications(self):
        self.helper_test_known_output_svs(self.test_objects_simple_dups[0], ['CACA'])
        self.helper_test_known_output_svs(self.test_objects_simple_dups[1], ['CACAT', 'CATAT'])
        self.helper_test_known_output_svs(self.test_objects_simple_dups[2], ['CC'])

    def test_simple_insertions(self):
        # test for INS with random insertion sequences -- going to check the post-mutation length and for the
        # correct placement of the original characters
        changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_simple_inss[0])
        hap_bools = [len(frag) == 7 and (frag[0] == 'C' or frag[-1] == 'A') for frag in
                     [changed_frag_1, changed_frag_2]]
        self.assertTrue(any(hap_bools))

    def test_simple_inversions(self):
        self.helper_test_known_output_svs(self.test_objects_simple_invs[0], ['TG'])
        self.helper_test_known_output_svs(self.test_objects_simple_invs[1], ['G'])

    def test_bidirectional_dispersion_events(self):
        # TRA
        # same output for forward and backward TRA
        self.helper_test_known_output_svs(self.test_dispersion_objects[0], ['TC'])
        # dDUP
        # different output for forward and backward dDUP
        self.helper_test_known_output_svs(self.test_dispersion_objects[1], ['CTC', 'TCT'])
        # INV_dDUP
        self.helper_test_known_output_svs(self.test_dispersion_objects[2], ['CTG', 'ACT'])
        # dDUP_iDEL
        self.helper_test_known_output_svs(self.test_dispersion_objects[3], ['CTC', 'GTG'])
        # INS_iDEL
        self.helper_test_known_output_svs(self.test_dispersion_objects[4], ['TC', 'GT'])

    def test_overlap_placement(self):
        # simple events
        for i in range(len(self.test_objects_overlap_simple)):
            config = self.test_objects_overlap_simple[i]
            config.initialize_files()
            curr_sim = SV_Simulator(config.ref, config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
            if i == 0:
                self.assertTrue('CTGTCGTA' in [changed_frag_1, changed_frag_2])
            elif i == 1:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 1)
                self.assertTrue('CCCC' in changed_frag_1 or 'CCCC' in changed_frag_2)
            elif i == 2:
                self.assertTrue('AA' not in changed_frag_1 or 'AA' not in changed_frag_2)
            elif i == 3:
                self.assertIsNone(curr_sim.svs[0].overlap_event)
            elif i == 4:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 0)
            elif i == 5:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 1)
            elif i == 6:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 0)
        # complex events
        for i in range(len(self.test_objects_overlap_cplx)):
            if i == 0:
                self.helper_test_known_output_svs(self.test_objects_overlap_cplx[i], ['CTGATGA', 'CGATGAT'])
            elif i == 1:
                # TRA [5,10) -> [10]; INV_dDUP [1,2) -> [0] or [3] - source ref: CTGATATGGAC
                changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_overlap_cplx[i])
                # need to account for the events being placed on opposite haplotypes, so will check for each separately
                # --> check for INV_dDUP in first four characters of output refs
                self.assertTrue(changed_frag_1[:4] in ['CTGA', 'ACTG'] or changed_frag_2[:4] in ['CTGA', 'ACTG'])
                # --> check for TRA in second half of refs
                self.assertTrue(changed_frag_1[-7:] in ['TCATGGA', 'ATGGATC'] or changed_frag_2[-7:] in ['TCATGGA', 'ATGGATC'])

    def test_divergence_events(self):
        # the divergence operator will mutate each base in an event interval with probability p
        # --> going to check for randomized placement of a divergence by checking that the output sequence
        # --> is not contained in the unedited reference (for event of length 5 and dummy reference: CTCCGTCGTA)
        changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_divergence_event[0])
        self.assertTrue(changed_frag_1 not in self.test_objects_divergence_event[0].ref or
                        changed_frag_2 not in self.test_objects_divergence_event[0].ref)

    def test_flanked_inversions(self):
        # tests for dupINVdup, delINVdel, etc.
        # --> all getting the reference CTCCGT
        targets = {'dupINVdup': 'CTACGGAGGT',
                   'delINVdel': 'GG',
                   'delINVdup': 'ACGGGT',
                   'dupINVdel': 'CTGGAG'}
        for i in range(4):
            self.helper_test_known_output_svs(self.test_objects_no_dis[i + 5], [targets[list(targets.keys())[i]]])

    def test_delINV_INVdel_events(self):
        # starter reference CTCCGT; each event of length 3
        targets = {'delINV': 'ACG',
                   'INVdel': 'GAG'}
        for i in range(2):
            self.helper_test_known_output_svs(self.test_objects_no_dis[i + 9], [targets[list(targets.keys())[i]]])

    def test_inverted_duplication_events(self):
        self.helper_test_known_output_svs(self.test_objects_no_dis[11], ['ACGACG'])

    def test_avoid_intervals(self):
        # test of calling choose_rand_pos and avoiding pre-specified intervals
        config_no_dis = self.test_objects_no_dis[0]
        config_no_dis.initialize_files()
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par)
        curr_sim.sim_settings['max_tries'] = 2000
        # define an interval in the simulation events dict that will be avoided in the pos choosing procedure
        curr_sim.event_ranges["Chromosome19"] = [(15, 83)]
        # the sv to be added is 15 bases long -- must go in spots 0-15
        curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)
        self.assertEqual(curr_sim.event_ranges, {'Chromosome19': [(15, 83), (0, 15)]})

    def test_load_avoid_intervals(self):
        # test of avoiding pre-specified intervals by reading them in from a vcf
        config_no_dis = self.test_objects_no_dis[4]
        config_no_dis.initialize_files()
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par)
        # --> example config is a chr19 DEL from pos 15-83 (same as above test)
        curr_sim.sim_settings['max_tries'] = 2000
        curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)
        self.assertEqual(curr_sim.event_ranges, {'chr19': [(15, 83), (0, 15)]})

    def test_chromosome_filtering(self):
        for i in range(len(self.test_objects_filter_chroms)):
            config_filter_chrom = self.test_objects_filter_chroms[i]
            config_filter_chrom.initialize_files()
            curr_sim = SV_Simulator(config_filter_chrom.ref, config_filter_chrom.par)
            self.assertEqual(set(curr_sim.len_dict.keys()), ({'chr21'} if i == 0 else set()))

    def test_req_space_filtering(self):
        config_req_sapce = self.test_objects_req_space[0]
        config_req_sapce.initialize_files()
        curr_sim = SV_Simulator(config_req_sapce.ref, config_req_sapce.par)
        with self.assertRaises(Exception):
            curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)


if __name__ == "__main__":
    unittest.main()
