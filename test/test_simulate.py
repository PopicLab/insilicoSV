import unittest
import sys
import os

from insilicosv.simulate import SV_Simulator
from pysam import FastaFile
import yaml
from insilicosv import utils
import numpy as np


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
        self.ref_file = "test/inputs/test.fna"
        self.par = "test/inputs/par.yaml"

        self.hap1 = "test/inputs/test1.fna"
        self.hap2 = "test/inputs/test2.fna"
        self.bed = "test/inputs/out.bed"

        self.test_overlap_bed = "test/inputs/example_overlap_events.bed"
        self.test_overlap_bed_2 = "test/inputs/example_overlap_events_2.bed"
        self.test_overlap_bed_4 = "test/inputs/example_overlap_events_4.bed"
        self.test_overlap_bed_5 = "test/inputs/example_overlap_events_5.bed"
        self.test_overlap_bed_6 = "test/inputs/example_overlap_events_6.bed"
        self.test_overlap_bed_7 = "test/inputs/example_overlap_events_7.bed"
        self.test_overlap_bed_8 = "test/inputs/example_overlap_events_8.bed"
        self.test_overlap_bed_9 = "test/inputs/example_overlap_events_9.bed"
        self.test_overlap_bed_10 = "test/inputs/example_overlap_events_10.bed"
        self.test_overlap_bed_11 = "test/inputs/example_overlap_events_11.bed"
        self.test_overlap_bed_12 = "test/inputs/example_overlap_events_12.bed"
        self.test_overlap_bed_13 = "test/inputs/example_overlap_events_13.bed"

        self.test_objects_no_dis = [TestObject([self.ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "delINVdup", "number": 1, "max_length": [5],
                                                    "min_length": [5]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "delINVdup", "number": 1, "max_length": [5], "min_length": [5]},
                                                   {"type": "delINVdel", "number": 1, "min_length": [5], "max_length": [5]},
                                                   {"type": "dupINVdup", "number": 1, "min_length": [5], "max_length": [5]}
                                                   ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGAGTCAGGGAGCAAAAAAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "dupINVdel", "number": 1, "max_length": [5], "min_length": [5]},
                                                   {"type": "delINV", "number": 1, "min_length": [5], "max_length": [5]},
                                                   {"type": "INVdel", "number": 1, "min_length": [5], "max_length": [5]}
                                                   ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "dup_INV", "number": 1, "max_length": [5], "min_length": [5]},
                                                   {"type": "INV_dup", "number": 1, "min_length": [5], "max_length": [5]}
                                                   ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "avoid_intervals": "test/inputs/example_avoid_interval.vcf",
                                                           "variant_sets": [{"type": "delINVdup", "number": 1,
                                                                    "max_length": [5], "min_length": [5]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # small ref for testing three-part events
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "dupINVdup",
                                                    "number": 1,
                                                    "min_length": [2, 2, 2],
                                                    "max_length": [2, 2, 2]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "delINVdel",
                                                    "number": 1,
                                                    "min_length": [2, 2, 2],
                                                    "max_length": [2, 2, 2]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "delINVdup",
                                                    "number": 1,
                                                    "min_length": [2, 2, 2],
                                                    "max_length": [2, 2, 2]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "dupINVdel",
                                                    "number": 1,
                                                    "min_length": [2, 2, 2],
                                                    "max_length": [2, 2, 2]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # objects for delINV and INVdel
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "delINV",
                                                    "number": 1,
                                                    "min_length": [3, 3],
                                                    "max_length": [3, 3]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "INVdel",
                                                    "number": 1,
                                                    "min_length": [3, 3],
                                                    "max_length": [3, 3]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # object for inverted duplication
                                    TestObject([self.ref_file, {"Chromosome19": "CGT"}],
                                               [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                   {"type": "INVdup",
                                                    "number": 1,
                                                    "min_length": [3],
                                                    "max_length": [3]}]}],
                                               self.hap1, self.hap2, self.bed)
                                    ]
        # test objects for bidirectional tests
        self.test_dispersion_objects = [TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                       {"type": "TRA",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                       {"type": "dDUP",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                       {"type": "INV_dDUP",
                                                        "number": 1,
                                                        "min_length": [1, 1],
                                                        "max_length": [1, 1]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTG"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                       {"type": "dDUP_iDEL",
                                                        "number": 1,
                                                        "min_length": [1, 1, 1],
                                                        "max_length": [1, 1, 1]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTG"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                       {"type": "INS_iDEL",
                                                        "number": 1,
                                                        "min_length": [1, 1, 1],
                                                        "max_length": [1, 1, 1]}]}],
                                                   self.hap1, self.hap2, self.bed)
                                        ]
        self.test_objects_ins = [TestObject([self.ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                            [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True}, "variant_sets": [
                                                {"type": "INS", "number": 1, "min_length": [5], "max_length": [5]},
                                                {"type": "delINV", "number": 1, "min_length": [5], "max_length": [5]},
                                                {"type": "INS", "number": 1, "min_length": [5], "max_length": [5]}]}],
                                            self.hap1, self.hap2, self.bed)]

        # --------- simple event test objects -----------
        self.test_objects_simple_dels = [TestObject([self.ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "DEL", "number": 1,
                                                                    "min_length": [13], "max_length": [13]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "DEL", "number": 1,
                                                                    "min_length": [14], "max_length": [14]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_dups = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "DUP", "number": 1,
                                                                    "min_length": [2], "max_length": [2]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "CAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "DUP", "number": 1,
                                                                    "min_length": [2], "max_length": [2]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "C"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "DUP", "number": 1,
                                                                    "min_length": [1], "max_length": [1]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_inss = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "INS", "number": 1,
                                                                    "min_length": [5], "max_length": [5]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_invs = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "INV", "number": 1,
                                                                    "min_length": [2], "max_length": [2]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "C"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                           "variant_sets": [{"type": "INV", "number": 1,
                                                                    "min_length": [1], "max_length": [1]}]}],
                                                    self.hap1, self.hap2, self.bed)]
        # ---------- test objects for overlap-aware event placement ------------
        self.test_objects_overlap_simple = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {"bed": self.test_overlap_bed},
                                                              "variant_sets": [{"type": "DEL", "number": 1,
                                                                       "min_length": [2], "max_length": [2],
                                                                       "num_overlap": 1}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {"bed": self.test_overlap_bed},
                                                              "variant_sets": [{"type": "DUP", "number": 1,
                                                                       "min_length": [2], "max_length": [2],
                                                                       "num_overlap": 1},
                                                                      {"type": "INV", "number": 1,
                                                                       "min_length": [2], "max_length": [2],
                                                                       "num_overlap": 1}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # combine two input files, filter all but one event by type
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {"bed": [self.test_overlap_bed, self.test_overlap_bed_2],
                                                                                 "allow_types": ["L1PA15"]},
                                                              "variant_sets": [{"type": "DEL", "number": 1,
                                                                       "min_length": [2], "max_length": [5],
                                                                       "num_overlap": 1}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # combine two input files, filter all by length
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {"bed": [self.test_overlap_bed, self.test_overlap_bed_2]},
                                                              "variant_sets": [{"type": "DEL", "number": 1,
                                                                       "min_length": [10], "max_length": [10],
                                                                       "num_overlap": 1}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # combine two input files, filter all by chromosome
                                            TestObject([self.ref_file, {"chr19": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [self.test_overlap_bed, self.test_overlap_bed_2]},
                                                              "variant_sets": [{"type": "DEL", "number": 1,
                                                                       "min_length": [10], "max_length": [10],
                                                                       "num_overlap": 1}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # type-specific num_overlap params
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [self.test_overlap_bed, self.test_overlap_bed_2],
                                                                  "allow_types": ["L1HS", "ALR/Alpha"]},
                                                              "variant_sets": [{"type": "DEL", "number": 4,
                                                                       "min_length": [1], "max_length": [10],
                                                                       "num_overlap": [2, 1]}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # type-specific num_overlap param > num available (ALR)
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {
                                                                  "bed": [self.test_overlap_bed, self.test_overlap_bed_2],
                                                                  "allow_types": ["L1HS", "ALR/Alpha"]},
                                                              "variant_sets": [{"type": "DEL", "number": 5,
                                                                       "min_length": [1], "max_length": [5],
                                                                       "num_overlap": [2, 3]}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            TestObject([self.ref_file, {"chr21": "CCTCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTATCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                              "overlap_events": {"bed": self.test_overlap_bed_11,
                                                                                 "allow_types": ['Alu', 'L1', 'L2', 'SVA', 'HERVK']},
                                                              "variant_sets": [{"type": "DEL", "number": 5,
                                                                       "min_length": [2], "max_length": [4],
                                                                       "num_overlap": [1, 1, 1, 1, 1]},
                                                                      {"type": "DEL", "number": 5,
                                                                       "min_length": [6], "max_length": [8],
                                                                       "num_overlap": [1, 1, 1, 1, 1]}
                                                                      ]}],
                                                       self.hap1, self.hap2, self.bed)
                                            ]
        self.test_objects_overlap_cplx = [TestObject([self.ref_file, {"chr21": "CTGAT"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                            "overlap_events": {"bed": [self.test_overlap_bed, self.test_overlap_bed_2],
                                                                               "allow_types": ["L1HS"]},
                                                            "variant_sets": [{"type": "dDUP", "number": 1,
                                                                     "min_length": [2, 1],
                                                                     "max_length": [2, 1],
                                                                     "num_overlap": 1}]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTGATATGGAC"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                            "overlap_events": {
                                                                "bed": [self.test_overlap_bed, self.test_overlap_bed_2],
                                                                "allow_types": ["L1HS", "AluSz6"]},
                                                            "variant_sets": [{"type": "TRA", "number": 1,
                                                                     "min_length": [4, 1],
                                                                     "max_length": [6, 1],
                                                                     "num_overlap": 1},
                                                                    {"type": "INV_dDUP", "number": 1,
                                                                     "min_length": [1, 1],
                                                                     "max_length": [1, 1],
                                                                     "num_overlap": 1}
                                                                    ]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                            "overlap_events": {"bed": self.test_overlap_bed_2},
                                                            "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                     "min_length": [3, 3, 3],
                                                                     "max_length": [3, 3, 3],
                                                                     "num_overlap": 1}]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                            "overlap_events": {"bed": self.test_overlap_bed_10},
                                                            "variant_sets": [{"type": "delINV", "number": 1,
                                                                     "min_length": [3, 3],
                                                                     "max_length": [3, 3],
                                                                     "num_overlap": 1},
                                                                    {"type": "INVdel", "number": 1,
                                                                     "min_length": [2, 2],
                                                                     "max_length": [2, 2],
                                                                     "num_overlap": 1}
                                                                    ]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                            "overlap_events": {"bed": self.test_overlap_bed_13},
                                                            "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                     "min_length": [3, 3, 3],
                                                                     "max_length": [3, 3, 3],
                                                                     "num_partial_overlap": 1}]}],
                                                     self.hap1, self.hap2, self.bed)
                                          ]
        self.test_objects_alu_mediated = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                             "fail_if_placement_issues": True},
                                                            "overlap_events": {"bed": self.test_overlap_bed_4},
                                                            "variant_sets": [{"type": "DEL", "number": 1,
                                                                     "min_length": [13], "max_length": [15],
                                                                     "num_alu_mediated": 1}]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTCCTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                             "fail_if_placement_issues": True},
                                                            "overlap_events": {"bed": [self.test_overlap_bed_4,
                                                                                       self.test_overlap_bed_5]},
                                                            "variant_sets": [{"type": "DEL", "number": 1,
                                                                     "min_length": [14], "max_length": [14],
                                                                     "num_alu_mediated": 1},
                                                                    {"type": "DEL", "number": 1,
                                                                     "min_length": [5], "max_length": [7],
                                                                     "num_alu_mediated": 1}
                                                                    ]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                             "fail_if_placement_issues": True},
                                                            "overlap_events": {"bed": [self.test_overlap_bed_4,
                                                                                       self.test_overlap_bed_5]},
                                                            "variant_sets": [{"type": "DEL", "number": 5,
                                                                     "min_length": [2], "max_length": [2],
                                                                     "num_alu_mediated": 5}]}],
                                                     self.hap1, self.hap2, self.bed)
                                          ]
        self.test_objects_known_elt_mix = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": True},
                                                             "overlap_events": {"bed": [self.test_overlap_bed, self.test_overlap_bed_4],
                                                                                "allow_types": "L1HS"},
                                                             "variant_sets": [{"type": "DEL", "number": 1,
                                                                      "min_length": [13], "max_length": [15],
                                                                      "num_alu_mediated": 1},
                                                                     {"type": "dDUP", "number": 1,
                                                                      "min_length": [2, 1],
                                                                      "max_length": [2, 1],
                                                                      "num_overlap": 1}]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": True},
                                                             "overlap_events": {
                                                                 "bed": [self.test_overlap_bed_6, self.test_overlap_bed_4]},
                                                             "variant_sets": [{"type": "DEL", "number": 1,
                                                                      "min_length": [13], "max_length": [15],
                                                                      "num_alu_mediated": 1},
                                                                     {"type": "dDUP", "number": 1,
                                                                      "min_length": [2, 1],
                                                                      "max_length": [2, 1],
                                                                      "num_overlap": 1}]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": False},
                                                             "overlap_events": {"bed": self.test_overlap_bed_4,
                                                                                "allow_types": "Alu"},
                                                             "variant_sets": [{"type": "DEL", "number": 1,
                                                                      "min_length": [4], "max_length": [4],
                                                                      "num_overlap": 1},
                                                                     {"type": "DEL", "number": 1,
                                                                      "min_length": [13], "max_length": [15],
                                                                      "num_alu_mediated": 1}
                                                                     ]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": False},
                                                             "overlap_events": {"bed": [self.test_overlap_bed_4,
                                                                                        self.test_overlap_bed_5],
                                                                                "allow_types": ["Alu", "L1", "MLT", "(GT)"]},
                                                             "variant_sets": [{"type": "DEL", "number": 3,
                                                                      "min_length": [2], "max_length": [8],
                                                                      "num_overlap": [0, 2, 1, 0]}
                                                                     ]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": False},
                                                             "overlap_events": {"bed": [self.test_overlap_bed_2,
                                                                                        self.test_overlap_bed_4,
                                                                                        self.test_overlap_bed_5],
                                                                                "allow_types": ["Alu", "ALR", "L1", "MLT", "(GT)"]},
                                                             "variant_sets": [{"type": "DEL", "number": 5,
                                                                      "min_length": [2], "max_length": [8],
                                                                      "num_overlap": [0, 1, 2, 0, 0],
                                                                      "num_partial_overlap": [0, 1, 1, 0, 0]}
                                                                     ]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "fail_if_placement_issues": False},
                                                             "overlap_events": {"bed": self.test_overlap_bed_9,
                                                                                "allow_types": ["Alu", "ALR", "L1", "MLT", "(GT)"]},
                                                             "variant_sets": [{"type": "DEL", "number": 6,
                                                                      "min_length": [2], "max_length": [4],
                                                                      "num_overlap": [0, 1, 2, 0, 0],
                                                                      "num_partial_overlap": [0, 1, 0, 0, 0],
                                                                      "num_alu_mediated": 1}
                                                                     ]}],
                                                      self.hap1, self.hap2, self.bed)
                                           ]
        self.test_objects_partial_overlap = [TestObject([self.ref_file, {"chr21": "CTCCGTAGTA"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                                "fail_if_placement_issues": True},
                                                               "overlap_events": {"bed": self.test_overlap_bed_12},
                                                               "variant_sets": [{"type": "DEL", "number": 1,
                                                                        "min_length": [2], "max_length": [2],
                                                                        "num_partial_overlap": 1}]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAG"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                                "fail_if_placement_issues": True},
                                                               "overlap_events": {"bed": self.test_overlap_bed_7},
                                                               "variant_sets": [{"type": "DEL", "number": 2,
                                                                        "min_length": [2], "max_length": [6],
                                                                        "num_partial_overlap": 1,
                                                                        "num_overlap": 1}]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAGGTCTAGC"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                                "fail_if_placement_issues": True},
                                                               "overlap_events": {"bed": self.test_overlap_bed_8},
                                                               "variant_sets": [{"type": "DEL", "number": 3,
                                                                        "min_length": [2], "max_length": [6],
                                                                        "num_partial_overlap": 1,
                                                                        "num_overlap": 1,
                                                                        "num_alu_mediated": 1}]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAGGTCTAGC"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                                "fail_if_placement_issues": True},
                                                               "overlap_events": {"bed": self.test_overlap_bed_8},
                                                               "variant_sets": [{"type": "DEL", "number": 2,
                                                                        "min_length": [2], "max_length": [6],
                                                                        "num_partial_overlap": 1,
                                                                        "num_alu_mediated": 1}]}],
                                                        self.hap1, self.hap2, self.bed)
                                             ]

        # ---------- test objects for divergence event ------------
        self.test_objects_divergence_event = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                                "variant_sets": [{"type": "DIVERGENCE", "number": 1,
                                                                         "min_length": [5], "max_length": [5]}]}],
                                                         self.hap1, self.hap2, self.bed),
                                              TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                                "variant_sets": [{"type": "DIVERGENCE", "number": 1, 'divergence_prob': 1,
                                                                         "min_length": [10], "max_length": [10]}]}],
                                                         self.hap1, self.hap2, self.bed),
                                              TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                                "variant_sets": [{"type": "DIVERGENCE", "number": 1, 'divergence_prob': 0.2,
                                                                         "min_length": [10], "max_length": [10]}]}],
                                                         self.hap1, self.hap2, self.bed)
                                              ]

        self.test_objects_filter_chroms = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                                                  "chr20": "CTCCGT"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "filter_small_chr": 10},
                                                             "variant_sets": [{"type": "DEL", "number": 1,
                                                                      "max_length": [3],
                                                                      "min_length": [3]}]}],
                                                      self.hap1, self.hap2, self.bed),
                                           TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                                                  "chr20": "CTCCGT"}],
                                                      [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True,
                                                                              "filter_small_chr": 50},
                                                             "variant_sets": [{"type": "DEL", "number": 1,
                                                                      "max_length": [3],
                                                                      "min_length": [3]}]}],
                                                      self.hap1, self.hap2, self.bed)]

        self.test_objects_req_space = [TestObject([self.ref_file, {"chr21": "CTCCGT"}],
                                                  [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                         "variant_sets": [{"type": "DEL", "number": 1,
                                                                  "max_length": [9],
                                                                  "min_length": [9]}]}],
                                                  self.hap1, self.hap2, self.bed)]
        # ----------- test objects for SNPs -----------
        self.test_objects_SNPs = [TestObject([self.ref_file, {"chr21": "C"}],
                                             [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                    "variant_sets": [{"type": "SNP", "number": 1}]}],
                                             self.hap1, self.hap2, self.bed),
                                  TestObject([self.ref_file, {"chr21": "CTGTTGACCG"}],
                                             [self.par, {"sim_settings": {"reference": self.ref_file, "prioritize_top": True},
                                                    "variant_sets": [{"type": "SNP", "number": 4}]}],
                                             self.hap1, self.hap2, self.bed)
                                  ]

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
        curr_sim = SV_Simulator(config.par)
        curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
        changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
        config.remove_test_files()
        if target_frags is not None:
            self.assertTrue(changed_frag_1 in target_frags or changed_frag_2 in target_frags)
        return changed_frag_1, changed_frag_2

    def test_snps(self):
        self.helper_test_known_output_svs(self.test_objects_SNPs[0], ['A', 'G', 'T'])
        changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_SNPs[1])
        self.assertTrue(changed_frag_1 not in self.test_objects_SNPs[1].ref or
                        changed_frag_2 not in self.test_objects_SNPs[1].ref)
        self.assertTrue(len(changed_frag_1) == len(changed_frag_2) == 10)

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

    def test_overlap_placement_simple(self):
        # simple events
        for i in range(len(self.test_objects_overlap_simple)):
            if i != 7:
                continue
            config = self.test_objects_overlap_simple[i]
            config.initialize_files()
            curr_sim = SV_Simulator(config.par)
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
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 2)
            elif i == 6:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 1)
            elif i == 7:
                self.assertEqual(len(curr_sim.overlap_events.overlap_events_dict.values()), 3)

    def test_overlap_placement_complex(self):
        # complex events
        for i in range(len(self.test_objects_overlap_cplx)):
            if i != 4:
                continue
            config = self.test_objects_overlap_cplx[i]
            config.initialize_files()
            curr_sim = SV_Simulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            # changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
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
            elif i == 2:
                for sv in curr_sim.svs:
                    self.assertIn([(ev.start, ev.end) for ev in sv.source_events],
                                  [[(10 + i * 3, 13 + i * 3), (13 + i * 3, 16 + i * 3), (16 + i * 3, 19 + i * 3)] for i in range(3)])
            elif i == 3:
                for sv in curr_sim.svs:
                    if sv.type == 'delINV':
                        self.assertIn([(ev.start, ev.end) for ev in sv.source_events],
                                      [[(15 + i * 3, 18 + i * 3), (18 + i * 3, 21 + i * 3)] for i in range(3)])
                    elif sv.type == 'INVdel':
                        self.assertIn([(ev.start, ev.end) for ev in sv.source_events],
                                      [[(8 + i * 2, 10 + i * 2), (10 + i * 2, 12 + i * 2)] for i in range(3)])
            elif i == 4:
                # partial overlap for flanked INV
                possible_shifts = [(10 + i * 3, 13 + i * 3, 16 + i * 3) for i in range(3)]
                correct_shift = False
                (frag_a, frag_b, frag_c) = tuple([ev.start for ev in curr_sim.svs[0].source_events])
                for (ivl_a, ivl_b, ivl_c) in possible_shifts:
                    if all(np.abs(frag - ivl) < 3 for (frag, ivl) in zip((ivl_a, ivl_b, ivl_c), (frag_a, frag_b, frag_c))):
                        correct_shift = True
                self.assertTrue(correct_shift)

    def test_partial_overlap_placement(self):
        for i in range(len(self.test_objects_partial_overlap)):
            config = self.test_objects_partial_overlap[i]
            config.initialize_files()
            curr_sim = SV_Simulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            # to simplify checking the diffnt possible outcomes of each test: sorting SVs by start position
            curr_sim.svs.sort(key=lambda x: x.start)
            changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
            if i == 0:
                # source: CTCCGTAGTA -> four possible valid del intervals, checking each
                self.assertTrue(any(possible_del not in changed_frag_1 for possible_del in ['CT', 'CC', 'CG', 'GT']) or
                                any(possible_del not in changed_frag_2 for possible_del in ['CT', 'CC', 'CG', 'GT']))
            if i == 1:
                # check that between the two SVs one is a perfect match and is other overlaps
                case_a = (curr_sim.svs[0].start, curr_sim.svs[0].end) == (2, 4) and utils.is_overlapping([(15, 20)], (curr_sim.svs[1].start, curr_sim.svs[1].end))
                case_b = (curr_sim.svs[1].start, curr_sim.svs[1].end) == (15, 20) and utils.is_overlapping([(2, 4)], (curr_sim.svs[0].start, curr_sim.svs[0].end))
                self.assertTrue(case_a or case_b)
            if i == 2:
                case_a = (curr_sim.svs[0].start, curr_sim.svs[0].end) == (2, 4) and utils.is_overlapping([(22, 25)], (curr_sim.svs[2].start, curr_sim.svs[2].end))
                case_b = (curr_sim.svs[2].start, curr_sim.svs[2].end) == (22, 25) and utils.is_overlapping([(2, 4)], (curr_sim.svs[0].start, curr_sim.svs[0].end))
                self.assertTrue(case_a or case_b)
                self.assertTrue((curr_sim.svs[1].start, curr_sim.svs[1].end) == (11, 15))
            if i == 3:
                sv_intervals = [(sv.start, sv.end) for sv in curr_sim.svs]
                self.assertTrue((11, 15) in sv_intervals)
                partial_ovl_sv = curr_sim.svs[(sv_intervals.index((11, 15)) + 1) % 2]
                self.assertTrue(utils.is_overlapping([(2, 4)], (partial_ovl_sv.start, partial_ovl_sv.end)) or
                                utils.is_overlapping([(22, 25)], (partial_ovl_sv.start, partial_ovl_sv.end)))

    def test_alu_mediated_placement(self):
        for i in range(len(self.test_objects_alu_mediated)):
            changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_alu_mediated[i])
            # single alu-mediated interval selection
            if i == 0:
                self.assertTrue('CTCCGTCTCCGTCGTACTAAGTCGTA' in [changed_frag_1, changed_frag_2])
            # bimodal SV config given: selection of two alu-mediated intervals
            if i == 1:
                # source: CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA
                # --> del1: GTACTAAGTCGTAC; del2: CGTCCT; --> want to check that both are missing from at least one hap
                self.assertFalse(all('GTACTAAGTCGTAC' in frag for frag in [changed_frag_1, changed_frag_2]))
                self.assertFalse(all('CGTCCT' in frag for frag in [changed_frag_1, changed_frag_2]))
            # i == 2: number of available intervals < number of alu-mediated SVs requested
            # --> no assert statement: will raise an error if any of the SVs are failed to be placed

    def test_mixed_known_element_placement(self):
        # handling for interacting with alu-mediated placement and overlap placement with other known elements
        for i in range(len(self.test_objects_known_elt_mix)):
            config = self.test_objects_known_elt_mix[i]
            config.initialize_files()
            curr_sim = SV_Simulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
            sv_intervals = [(sv.start, sv.end) for sv in curr_sim.svs]
            if i == 0:
                # source: CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA
                # --> del: GTACTAAGTCGTAC; dDUP: CCGCC [forward] or CCTCC [backward]; --> want to check that the del
                # missing from at least one hap and that the dDUP is present in at least one hap
                self.assertFalse(all('GTACTAAGTCGTAC' in frag for frag in [changed_frag_1, changed_frag_2]))
                self.assertTrue(any('CCGCC' in frag for frag in [changed_frag_1, changed_frag_2]) or
                                any('CCTCC' in frag for frag in [changed_frag_1, changed_frag_2]))
            if i == 1:
                # --> del: same as above; dDUP: source either at (2,4) or (25,27)
                self.assertFalse(all('GTACTAAGTCGTAC' in frag for frag in [changed_frag_1, changed_frag_2]))
                self.assertTrue(any('CCGCC' in frag for frag in [changed_frag_1, changed_frag_2]) or
                                any('CCTCC' in frag for frag in [changed_frag_1, changed_frag_2]) or
                                any('TCGTC' in frag for frag in [changed_frag_1, changed_frag_2]))
            if i == 2:
                # test case for valid alu-pair selection taking precedence over single element overlap selection
                # --> the single overlap element should fail to find the (18,22) Alu because it should be claimed by
                # --> the alu-mediated DEL that will need that for the right breakpoint
                self.assertIn((6, 20), sv_intervals)
                self.assertNotIn((18, 22), sv_intervals)
            if i == 3:
                # test case of specifying counts of zero for certain element types
                self.assertIn((10, 18), sv_intervals)
                self.assertIn((1, 4), sv_intervals)
                self.assertIn((23, 25), sv_intervals)
            if i == 4:
                # test case of specifying counts of zero for certain element types: split across full and partial ovlps
                # --> check for the right number of full overlaps of the specified types
                self.assertTrue(any(alr in sv_intervals for alr in [(10, 12), (16, 19)]))
                self.assertTrue([l1 in sv_intervals for l1 in [(13, 15), (10, 18), (1, 4)]].count(True) in [2, 3])
                # --> and the right number of partial overlaps
                self.assertTrue(any(utils.is_overlapping(sv_intervals, alu, strictly_partial=True) for alu in [(10, 12), (16, 19)]) or
                                [alu in sv_intervals for alu in [(10, 12), (16, 19)]].count(True) == 2)
                self.assertTrue(any(utils.is_overlapping(sv_intervals, l1, strictly_partial=True) for l1 in [(13, 15), (10, 18), (1, 4)]) or
                                [l1 in sv_intervals for l1 in [(13, 15), (10, 18), (1, 4)]].count(True) == 3)
            if i == 5:
                # --> alu-mediated interval has a single option
                self.assertIn((6, 8), sv_intervals)
                # --> both L1s in fixed locations
                self.assertIn((1, 3), sv_intervals)
                self.assertIn((25, 27), sv_intervals)
                self.assertTrue(any(alr in sv_intervals for alr in [(10, 12), (16, 19)]))
                # --> and the right number of partial overlaps
                self.assertTrue(any(utils.is_overlapping(sv_intervals, alu, strictly_partial=True) for alu in [(10, 12), (16, 19)]) or
                                [alu in sv_intervals for alu in [(10, 12), (16, 19)]].count(True) == 2)

    def test_divergence_events(self):
        # the divergence operator will mutate each base in an event interval with probability p
        # --> going to check for randomized placement of a divergence by checking that the output sequence
        # --> is not contained in the unedited reference (for event of length 5 and dummy reference: CTCCGTCGTA)
        for i in range(len(self.test_objects_divergence_event)):
            changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_divergence_event[i])
            self.assertTrue(changed_frag_1 not in self.test_objects_divergence_event[i].ref or
                            changed_frag_2 not in self.test_objects_divergence_event[i].ref)
            self.assertTrue(len(changed_frag_1) == len(changed_frag_2) == 10)


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
        curr_sim = SV_Simulator(config_no_dis.par)
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
        curr_sim = SV_Simulator(config_no_dis.par)
        # --> example config is a chr19 DEL from pos 15-83 (same as above test)
        curr_sim.sim_settings['max_tries'] = 2000
        curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)
        self.assertEqual(curr_sim.event_ranges, {'chr19': [(15, 83), (0, 15)]})

    def test_chromosome_filtering(self):
        for i in range(len(self.test_objects_filter_chroms)):
            config_filter_chrom = self.test_objects_filter_chroms[i]
            config_filter_chrom.initialize_files()
            curr_sim = SV_Simulator(config_filter_chrom.par)
            self.assertEqual(set(curr_sim.len_dict.keys()), ({'chr21'} if i == 0 else set()))

    def test_req_space_filtering(self):
        config_req_sapce = self.test_objects_req_space[0]
        config_req_sapce.initialize_files()
        curr_sim = SV_Simulator(config_req_sapce.par)
        with self.assertRaises(Exception):
            curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)


if __name__ == "__main__":
    unittest.main()
