import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pytest
import yaml

from insilicosv.simulate import SVSimulator
from pysam import FastaFile
from intervaltree import Interval
from insilicosv import utils
from insilicosv.utils import Region, as_list, reverse_complement
from insilicosv.sv_defs import (
    Transform, TransformType,
    Operation, SV)

def is_overlapping(event_ranges, addition, called_from_helper=False, strictly_partial=False):
    # addition: tuple (start, end)
    # event_ranges: list containing tuples
    # checks if addition overlaps with any of the events already stored
    # 'strictly_partial' a toggle to return True if the interval overlaps a stored event, without being equal to it
    for event in event_ranges:
        if event[1] > addition[0] and event[0] < addition[1]:
            if called_from_helper:
                raise Exception("Overlap between {} and {}".format(event[0:2], addition[0:2]))
            else:
                return True if not strictly_partial else event != addition

    return False

def get_span(sv: SV) -> Region:
    assert sv.is_placed()

    regions = sorted(sv.get_regions())
    # TODO: update for inter-chromosomal operations
    chrom = regions[0].chrom
    assert all(region.chrom == chrom for region in regions)
    return Region(chrom=chrom, start=regions[0].start, end=regions[-1].end)

class TestObject():

    __test__ = False

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
        with FastaFile(self.hap1) as fasta_out_1, FastaFile(self.hap2) as fasta_out_2:
            def get_frag(fasta):
                result = [fasta.fetch(ref) for ref in fasta.references]
                return result[0] if len(result) == 1 else tuple(result)
            if return_haps == 'hap1':
                return get_frag(fasta_out_1)
            elif return_haps == 'hap2':
                return get_frag(fasta_out_2)
            else:
                return (get_frag(fasta_out_1), get_frag(fasta_out_2))

class TestSVSimulator(unittest.TestCase):
    def setUp(self):

        self.test_dir = tempfile.mkdtemp()
        self.ref_file = f"{self.test_dir}/test.fna"
        self.par = f"{self.test_dir}/par.yaml"

        self.hap1 = f"{self.test_dir}/test1.fna"
        self.hap2 = f"{self.test_dir}/test2.fna"
        self.bed = f"{self.test_dir}/out.bed"

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

        self.test_exclude_bed = "test/inputs/exclude.bed"

        self.test_objects_no_dis = [TestObject([self.ref_file, {
            "chr21": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file,
                                                                  "max_tries": 2000},
                                                 "blacklist_regions": [self.test_exclude_bed],
                                                 "variant_sets": [
                                                     {"type": "delINVdup", "number": 1, 
                                                      "length_ranges": [[5, 5], [5, 5], [5, 5]],
                                                      "blacklist_region_type": "AluSz6"}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {
                                        "Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "delINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]},
                                                     {"type": "delINVdel", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]},
                                                     {"type": "dupINVdup", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]}
                                                 ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {
                                        "Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGAGTCAGGGAGCAAAAAAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "dupINVdel", "number": 1, "length_ranges": [[5, 5], [5, 5], [5, 5]]},
                                                     {"type": "delINV", "number": 1, "length_ranges": [[5, 5], [5, 5]]},
                                                     {"type": "INVdel", "number": 1, "length_ranges": [[5, 5], [5, 5]]}
                                                 ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {
                                        "Chromosome19": "ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "fldup_INV", "number": 1, "length_ranges": [[5, 5]]},
                                                     {"type": "INV_fldup", "number": 1, "length_ranges": [[5, 5]]}
                                                 ]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {
                                        "chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "blacklist_regions": "test/inputs/example_avoid_interval.vcf",
                                                 "variant_sets": [{"type": "delINVdup", "number": 1,
                                                                   "length_ranges": [[5, 5], [5, 5], [5, 5]],
                                                                   "blacklist_region_type": "all"}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # small ref for testing three-part events
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "dupINVdup",
                                                      "number": 1,
                                                      "length_ranges": [[2, 2], [2, 2], [2, 2]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "delINVdel",
                                                      "number": 1,
                                                      "length_ranges": [[2, 2], [2, 2], [2, 2]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "delINVdup",
                                                      "number": 1,
                                                      "length_ranges": [[2, 2], [2, 2], [2, 2]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "dupINVdel",
                                                      "number": 1,
                                                      "length_ranges": [[2, 2], [2, 2], [2, 2]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # objects for delINV and INVdel
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "delINV",
                                                      "number": 1,
                                                      "length_ranges": [[3, 3], [3, 3]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {"Chromosome19": "CTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "INVdel",
                                                      "number": 1,
                                                      "length_ranges": [[3, 3], [3, 3]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    # object for inverted duplication
                                    TestObject([self.ref_file, {"Chromosome19": "CGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "variant_sets": [
                                                     {"type": "INV_DUP3",
                                                      "number": 1,
                                                      "length_ranges": [[3, 3]]}]}],
                                               self.hap1, self.hap2, self.bed),
                                    TestObject([self.ref_file, {
                                        "chr19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                               [self.par,
                                                {"sim_settings": {"reference": self.ref_file},
                                                 "blacklist_regions": "test/inputs/example_avoid_interval.bed",
                                                 "variant_sets": [{"type": "delINVdup", "number": 1,
                                                                   "length_ranges": [[5, 5], [5, 5], [5, 5]],
                                                                   "blacklist_region_type": "all"}]}],
                                               self.hap1, self.hap2, self.bed)
                                    ]
        # test objects for bidirectional tests
        self.test_dispersion_objects = [TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file},
                                                               "variant_sets": [
                                                                   {"type": "TRA_NONRECIPROCAL",
                                                                    "number": 1,
                                                                    "length_ranges": [[1, 1], [1, 1]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file},
                                                               "variant_sets": [
                                                                   {"type": "dDUP",
                                                                    "number": 1,
                                                                    "length_ranges": [[1, 1], [1, 1]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CT"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file },
                                                               "variant_sets": [
                                                                   {"type": "INV_dDUP",
                                                                    "number": 1,
                                                                    "length_ranges": [[1, 1], [1, 1]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTG"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file },
                                                               "variant_sets": [
                                                                   {"type": "dDUP_iDEL",
                                                                    "number": 1,
                                                                    "length_ranges": [[1, 1], [1, 1], [1, 1]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTG"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file},
                                                               "variant_sets": [
                                                                   {"type": "INS_iDEL",
                                                                    "number": 1,
                                                                    "length_ranges": [[1, 1], [1, 1], [1, 1]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file}, "variant_sets": [
                                                       {"type": "TRA_RECIPROCAL",
                                                        "number": 1,
                                                        "length_ranges": [[1, 1], [1, 1], [3, 3]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file}, "variant_sets": [
                                                       {"type": "TRA_RECIPROCAL",
                                                        "number": 1,
                                                        "length_ranges": [[1, 1], [1, 1], [3, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file}, "variant_sets": [
                                                       {"type": "TRA_RECIPROCAL",
                                                        "number": 1,
                                                        "length_ranges": [[2, 2], [2, 1], [3, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file}, "variant_sets": [
                                                       {"type": "TRA_RECIPROCAL",
                                                        "number": 1,
                                                        "length_ranges": [[2, 2], [2, 2], [1, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file},
                                                               "variant_sets": [
                                                                   {"type": "TRA_RECIPROCAL",
                                                                    "number": 1,
                                                                    "length_ranges": [[3, 3], [2, 2], [None, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"Chromosome19": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file},
                                                               "variant_sets": [
                                                                   {"type": "TRA_RECIPROCAL",
                                                                    "number": 1,
                                                                    "length_ranges": [[2, 2], [2, 2], [None, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        TestObject([self.ref_file, {"ChromA": "CTTTA"}],
                                                   [self.par, {"sim_settings": {"reference": self.ref_file,
                                                                                "homozygous_only": True},
                                                               "variant_sets": [
                                                                   {"type": "TRA_RECIPROCAL",
                                                                    "number": 1,
                                                                    "length_ranges": [[2, 2], [2, 2], [None, None]]}]}],
                                                   self.hap1, self.hap2, self.bed),
                                        ]
        self.test_objects_ins = [TestObject([self.ref_file, {
            "Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}],
                                            [self.par,
                                             {"sim_settings": {"reference": self.ref_file},
                                              "variant_sets": [
                                                  {"type": "INS", "number": 1, "length_ranges": [[5, 5]]},
                                                  {"type": "delINV", "number": 1, "length_ranges": [[5, 5]]},
                                                  {"type": "INS", "number": 1, "length_ranges": [[5, 5]]}]}],
                                            self.hap1, self.hap2, self.bed)]

        # --------- simple event test objects -----------
        self.test_objects_simple_dels = [TestObject([self.ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "DEL", "number": 1,
                                                                                  "length_ranges": [[13, 13]]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "CACTATCTCTCCGAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "DEL", "number": 1,
                                                                                  "length_ranges": [[14, 14]]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_dups = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "DUP", "number": 1,
                                                                                  "length_ranges": [[2, 2]]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "CAT"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "DUP", "number": 1,
                                                                                  "length_ranges": [[2, 2]]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "C"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "DUP", "number": 1,
                                                                                  "length_ranges": [[1, 1]]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_inss = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "INS", "number": 1,
                                                                                  "length_ranges": [[5, 5]]}]}],
                                                    self.hap1, self.hap2, self.bed)]

        self.test_objects_simple_invs = [TestObject([self.ref_file, {"Chromosome19": "CA"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "INV", "number": 1,
                                                                                  "length_ranges": [[2, 2]]}]}],
                                                    self.hap1, self.hap2, self.bed),
                                         TestObject([self.ref_file, {"Chromosome19": "C"}],
                                                    [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                "variant_sets": [{"type": "INV", "number": 1,
                                                                                  "length_ranges": [[1, 1]]}]}],
                                                    self.hap1, self.hap2, self.bed)]
        # ---------- test objects for overlap-aware event placement ------------
        self.test_objects_overlap_simple = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                   "overlap_regions": [self.test_overlap_bed],
                                                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                                                     "length_ranges": [[None, None]],
                                                                                     "overlap_region_length_range": [2, 2],
                                                                                     "overlap_mode": "exact"}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                   "overlap_regions": [self.test_overlap_bed],
                                                                   "variant_sets": [{"type": "DUP", "number": 1,
                                                                                     "length_ranges": [[None, None]],
                                                                                     "overlap_region_length_range": [2, 2],
                                                                                     "overlap_mode": "exact"}
                                                                                    ]}],
                                                       self.hap1, self.hap2, self.bed),
                                            # combine two input files, filter all but one event by type
                                            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                                                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                   "overlap_regions": [self.test_overlap_bed,
                                                                                       self.test_overlap_bed_2],
                                                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                                                     "overlap_region_length_range": [2, 5],
                                                                                     "length_ranges": [[None, None]],
                                                                                     "overlap_mode": "exact",
                                                                                     "overlap_region_type": ["L1PA15"]}]}],
                                                       self.hap1, self.hap2, self.bed),
                                            ]
        self.test_objects_overlap_cplx = [TestObject([self.ref_file, {"chr21": "CTGAT"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "dDUP", "number": 1,
                                                                                   "length_ranges": [[None, None], [1, 1]],
                                                                                   "overlap_region_length_range": [2, 2],
                                                                                   "overlap_mode": "exact",
                                                                                   "overlap_region_type": ["L1HS"],
                                                                                   "source": "(A)_",
                                                                                   "target": "A_A^"}]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CTGATATGGAC"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "TRA_NONRECIPROCAL", "number": 1,
                                                                                   "length_ranges": [[None, None], [1, 1]],
                                                                                   "overlap_mode": "exact",
                                                                                   "overlap_region_length_range": [4, 6],
                                                                                   "overlap_region_type": ["L1HS"],
                                                                                   "source": "(A)_"},
                                                                                  {"type": "INV_dDUP", "number": 1,
                                                                                   "length_ranges": [[None, None], [1, 1]],
                                                                                   "overlap_mode": "exact",
                                                                                   "overlap_region_length_range": [1, 1],
                                                                                   "overlap_region_type": ["AluSz6"],
                                                                                   "source": "(A)_"
                                                                                   }]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file},
                                                "overlap_regions": [self.test_overlap_bed_2],
                                                "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                  "length_ranges": [[None, None], [3, 3], [3, 3]],
                                                                  "overlap_mode": "exact",
                                                                  "overlap_region_length_range": [3, 3],
                                                                  "source": "(A)BC"}]}],
                                              self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file},
                                                "overlap_regions": [self.test_overlap_bed_10],
                                                "variant_sets": [{"type": "delINV", "number": 1,
                                                                  "length_ranges": [[None, None], [3, 3]],
                                                                  "overlap_mode": "exact",
                                                                  "overlap_region_length_range": [3, 3],
                                                                  "overlap_region_type": ["ALR"],
                                                                  "source": "(A)B"},
                                                                 {"type": "INVdel", "number": 1,
                                                                  "length_ranges": [[None, None], [2, 2]],
                                                                  "overlap_region_length_range": [2, 2],
                                                                  "overlap_mode": "exact",
                                                                  "source": "(A)B"}
                                                                 ]}],
                                              self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file},
                                                "overlap_regions": [self.test_overlap_bed_13],
                                                "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                  "length_ranges": [[3, 3], [3, 3], [3, 3]],
                                                                  "overlap_mode": "partial",
                                                                  "source": "(A)BC"}]}],
                                              self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "dDUP", "number": 1,
                                                                                   "length_ranges": [[2, 2], [1, 1]],
                                                                                   "overlap_region_type": ["L1PA15"],
                                                                                   "overlap_mode": "contained",
                                                                                   "target": ["A", "_", "(", "A^", ")"]
                                                                                   }]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "INV_dDUP", "number": 1,
                                                                                   "length_ranges": [[2, 2], [1, 1]],
                                                                                   "overlap_region_type": ["L1PA15"],
                                                                                   "overlap_mode": "contained",
                                                                                   "target": ["A", "_", "(", "a^", ")"]
                                                                                   }]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CCTGATCTGATCTGATCTGATCTGATTGAT"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "TRA_NONRECIPROCAL", "number": 1,
                                                                                   "length_ranges": [[2, 2], [1, 1]],
                                                                                   "overlap_region_type": ["L1PA15"],
                                                                                   "overlap_mode": "contained",
                                                                                   "target": ["_", "(", "A^", ")"]
                                                                                   }]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject([self.ref_file, {"chr21": "CCTGATATGGACCTGATATGGACTGATATGGAC"}],
                                                     [self.par, {"sim_settings": {"reference": self.ref_file,
                                                                                  "max_tries": 1000},
                                                                 "overlap_regions": [self.test_overlap_bed,
                                                                                     self.test_overlap_bed_2],
                                                                 "variant_sets": [{"type": "INV_dDUP", "number": 1,
                                                                                   "length_ranges": [[2, 2], [3, 3]],
                                                                                   "overlap_region_type": ["ALR"],
                                                                                   "overlap_mode": "contained",
                                                                                   "target": ["A", "_", "(", "a^", ")"]},
                                                                                  {"type": "TRA_NONRECIPROCAL", "number": 1,
                                                                                   "length_ranges": [[4, 6], [1, 1]]}
                                                                                  ]}],
                                                     self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file,
                                                                 "max_tries": 1000},
                                                "overlap_regions": [self.test_overlap_bed_11],
                                                "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                  "length_ranges": [[2, 2], [3, 3], [2, 2]],
                                                                  "overlap_mode": "contained",
                                                                  "overlap_region_type": ["L1PA15"]},
                                                                 {"type": "DEL", "number": 1,
                                                                  "length_ranges": [[2, 6]]}
                                                                 ]}],
                                              self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file,
                                                                 "max_tries": 1000},
                                                "overlap_regions": [self.test_overlap_bed_11],
                                                "variant_sets": [{"type": "delINV", "number": 1,
                                                                  "length_ranges": [[3, 3], [4, 4]],
                                                                  "overlap_mode": "contained",
                                                                  "overlap_region_type": ["L1PA15"]}]}],
                                              self.hap1, self.hap2, self.bed),
                                          TestObject(
                                              [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                                              [self.par,
                                               {"sim_settings": {"reference": self.ref_file},
                                                "overlap_regions": [self.test_overlap_bed_11],
                                                "variant_sets": [{"type": "delINVdel", "number": 1,
                                                                  "length_ranges": [[2, 6], [3, 5], [2, 9]],
                                                                  "overlap_mode": "partial",
                                                                  "overlap_region_type": ["L1PA15"]}]}],
                                              self.hap1, self.hap2, self.bed),
                                          ]
        self.test_objects_frag_level_overlap = [
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[None, None], [1, 1]],
                                                     "overlap_region_length_range": [2, 2],
                                                     "overlap_mode": "exact",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "INV_dDUP", "number": 1,
                                                     "length_ranges": [[None, None], [1, 1]],
                                                     "overlap_region_length_range": [2, 2],
                                                     "overlap_mode": "exact",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "TRA_NONRECIPROCAL", "number": 1,
                                                     "length_ranges": [[None, None], [1, 1]],
                                                     "overlap_region_length_range": [2, 2],
                                                     "overlap_mode": "exact",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATTGATTGATGAGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_2],
                                   "variant_sets": [{"type": "dDUP_iDEL", "number": 1,
                                                     "length_ranges": [[None, None], [2, 2], [2, 2]],
                                                     "overlap_region_length_range": [2, 2],
                                                     "overlap_mode": "exact",
                                                     "source": ["(", "A", ")", "_", "B"],
                                                     "overlap_region_type": ["L1"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATTGATTGATGAA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_2],
                                   "variant_sets": [{"type": "dDUP_iDEL", "number": 1,
                                                     "length_ranges": [[2, 2], [None, None], [2, 2]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_length_range": [2, 2],
                                                     "source": ["A", "_", "(", "B", ")"],
                                                     "overlap_region_type": ["L1"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "partial",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "INV_dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "partial",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "TRA_NONRECIPROCAL", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "partial",
                                                     "source": ["(", "A", ")", "_"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "target": ["A", "_", "(", ")", "A^"],
                                                     "overlap_mode": "contained",
                                                     "overlap_region_type": ["L1HS"]
                                                     }]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "INV_dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "contained",
                                                     "target": ["A", "_", "(", "a^", ")"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "TRA_NONRECIPROCAL", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "target": ["_", "(", "A^", ")"],
                                                     "overlap_mode": "contained",
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "delINVdel", "number": 1,
                                                     "length_ranges": [[2, 2], [None, None], [2, 2]],
                                                     "overlap_region_length_range": [4, 6],
                                                     "overlap_mode": "exact",
                                                     "source": ["A", "(", "B", ")", "C"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "delINVdel", "number": 1,
                                                     "length_ranges": [[None, None], [2, 2], [2, 2]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_length_range": [4, 6],
                                                     "source": ["(", "A", ")", "B", "C"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "delINVdel", "number": 1,
                                                     "length_ranges": [[2, 2], [2, 2], [None, None]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_length_range": [4, 6],
                                                     "source": ["A", "B", "(", "C", ")"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "delINV", "number": 1,
                                                     "length_ranges": [[None, None], [2, 2]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_length_range": [4, 6],
                                                     "source": ["(", "A", ")", "B"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "INVdel", "number": 1,
                                                     "length_ranges": [[2, 2], [None, None]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_length_range": [4, 6],
                                                     "source": ["A", "(", "B", ")"],
                                                     "overlap_region_type": ["L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {
                "chr21": utils.generate_seq(100)}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [
                                       {"type": "Custom", "source": ["A", "(", "B", ")", "_", "C", "_", "D"],
                                        "target": ["b", "b^", "_", "A^", "E^", "c", "_", "E^", "D", "C^"], "number": 1,
                                        "length_ranges": [[2, 3], [None, None], [7, 10], [8, 10], [10, 15], [10, 15], [15, 20]],
                                        "overlap_mode": "exact",
                                        "overlap_region_length_range": [4, 10],
                                        "overlap_region_type": ["L1HS"]},
                                   ]}],
                       self.hap1, self.hap2, self.bed),
            ]

        self.test_objects_frag_level_overlap_unbounded_disps = [
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed,
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, None]],
                                                     "overlap_mode": "exact",
                                                     # want exact overlap with source, and target landing in the other
                                                     "overlap_region_type": ["L1HS", "L1HS"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed_2,
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 5], [1, None]],
                                                     "overlap_mode": "exact",
                                                     # want exact overlap with source, and target landing in the other
                                                     "overlap_region_type": ["L1", "Alu"]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTTGATGATGATTGATGATA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed_2,
                                   "variant_sets": [{"type": "dDUP", "number": 1,
                                                     "length_ranges": [[1, 5], [1, None]],
                                                     "overlap_mode": "exact",
                                                     # want exact overlap with source, and target landing in the other
                                                     "overlap_region_type": ["Alu", "L1"]}]}],
                       self.hap1, self.hap2, self.bed),
            ]
        self.test_objects_alu_mediated = [
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed_4,
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[13, 15]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTCCTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_4,
                                                       self.test_overlap_bed_5],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[14, 14]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"},
                                                    {"type": "DEL", "number": 1,
                                                     "length_ranges": [[5, 7]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"}
                                                    ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_4,
                                                       self.test_overlap_bed_5],
                                   "variant_sets": [{"type": "DEL", "number": 5,
                                                     "length_ranges": [[2, 2]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"}]}],
                       self.hap1, self.hap2, self.bed)
            ]
        self.test_objects_known_elt_mix = [
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed, self.test_overlap_bed_4],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[13, 15]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"},
                                                    {"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": "L1HS",
                                                     "overlap_component": "source"}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_6, self.test_overlap_bed_4],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[13, 15]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"},
                                                    {"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_mode": "exact",
                                                     "overlap_component": "source"}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_4,
                                                       self.test_overlap_bed_5],
                                   "variant_sets": [{"type": "DEL", "number": 2,
                                                     "length_ranges": [[2, 8]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": "L1"},
                                                    {"type": "DEL", "number": 1,
                                                     "length_ranges": [[2, 8]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": "MLT"}
                                                    ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed_9,
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[2, 4]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": "ALR"},
                                                    {"type": "DEL", "number": 2,
                                                     "length_ranges": [[2, 4]],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": "L1"},
                                                    {"type": "DEL", "number": 1,
                                                     "length_ranges": [[2, 4]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"},
                                                    {"type": "DEL", "number": 1,
                                                     "length_ranges": [[2, 4]]}
                                                    ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_2, self.test_overlap_bed_4],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[1, 5]],
                                                     "overlap_mode": "flanked",
                                                     "overlap_region_type": "Alu"},
                                                    {"type": "dDUP", "number": 1,
                                                     "length_ranges": [[2, 2], [1, 1]],
                                                     "overlap_region_type": "L1PA15",
                                                     "overlap_component": "target"}]}],
                       self.hap1, self.hap2, self.bed)
            ]
        self.test_objects_partial_overlap = [TestObject([self.ref_file, {"chr21": "CTCCGTAGTA"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                    "overlap_regions": [self.test_overlap_bed_12],
                                                                    "variant_sets": [{"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 2]],
                                                                                      "overlap_mode": "partial"}]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAG"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                    "overlap_regions": [self.test_overlap_bed_7],
                                                                    "variant_sets": [{"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "partial"},
                                                                                     {"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[None, None]],
                                                                                      "overlap_region_length_range": [2, 6],
                                                                                      "overlap_mode": "exact"}
                                                                                     ]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAGGTCTAGC"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                    "overlap_regions": [self.test_overlap_bed_8],
                                                                    "variant_sets": [{"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "partial"},
                                                                                     {"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "exact"},
                                                                                     {"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "flanked",
                                                                                      "overlap_region_type": "Alu"}
                                                                                     ]}],
                                                        self.hap1, self.hap2, self.bed),
                                             TestObject([self.ref_file, {"chr21": "CTCCGTAGTAAGTCAGGTGAGGCAGGTCTAGC"}],
                                                        [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                    "overlap_regions": [self.test_overlap_bed_8],
                                                                    "variant_sets": [{"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "partial"},
                                                                                     {"type": "DEL", "number": 1,
                                                                                      "length_ranges": [[2, 6]],
                                                                                      "overlap_mode": "flanked",
                                                                                      "overlap_region_type": "Alu"}
                                                                                     ]}],
                                                        self.hap1, self.hap2, self.bed)
                                             ]

        # ---------- test objects for divergence event ------------
        self.test_objects_divergence_event = [TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                     "variant_sets": [
                                                                         {"type": "DIVERGENCE", "number": 1,
                                                                          'divergence_prob': 0.5,
                                                                          "length_ranges": [[5, 5]]}]}],
                                                         self.hap1, self.hap2, self.bed),
                                              TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                     "variant_sets": [
                                                                         {"type": "DIVERGENCE", "number": 1,
                                                                          'divergence_prob': 1,
                                                                          "length_ranges": [[10, 10]]}]}],
                                                         self.hap1, self.hap2, self.bed),
                                              TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                                                         [self.par, {"sim_settings": {"reference": self.ref_file},
                                                                     "variant_sets": [
                                                                         {"type": "DIVERGENCE", "number": 1,
                                                                          'divergence_prob': 0.2,
                                                                          "length_ranges": [[10, 10]]}]}],
                                                         self.hap1, self.hap2, self.bed)
                                              ]

        self.test_objects_filter_chroms = [
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                        "chr20": "CTCCGT"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file,
                                                    "filter_small_chr": 10},
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[3, 3]]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA",
                                        "chr20": "CTCCGT"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file,
                                                    "filter_small_chr": 50},
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[3, 3]]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr1": "CTCCGT", "chr2": "CTCCGT", "chr3": "CTCCGT", "chrM": "C"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file,
                                                    "filter_small_chr": 4},
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "length_ranges": [[3, 3]]}]}],
                        self.hap1, self.hap2, self.bed)]

        self.test_objects_req_space = [TestObject([self.ref_file, {"chr21": "CTCCGT"}],
                                                  [self.par, {"sim_settings": {"reference": self.ref_file},
                                                              "variant_sets": [{"type": "DEL", "number": 1,
                                                                                "length_ranges": [[9, 9]]}]}],
                                                  self.hap1, self.hap2, self.bed)]
        # ----------- test objects for SNPs -----------
        self.test_objects_SNPs = [TestObject([self.ref_file, {"chr21": "C"}],
                                             [self.par,
                                              {"sim_settings": {"reference": self.ref_file},
                                               "variant_sets": [{"type": "SNP", "number": 1}]}],
                                             self.hap1, self.hap2, self.bed),
                                  TestObject([self.ref_file, {"chr21": "CTGTTGACCG"}],
                                             [self.par,
                                              {"sim_settings": {"reference": self.ref_file},
                                               "variant_sets": [{"type": "SNP", "number": 4}]}],
                                             self.hap1, self.hap2, self.bed)
                                  ]
        # ------------- test objects for blacklist regions
        self.test_blacklist_regions = [
            TestObject([self.ref_file, {
                "chr21": "ACTAATCTCTTCTCTCTTCTCTCTCCGT"}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file},
                         "blacklist_regions": "test/inputs/example_avoid_intervals_2.bed",
                         "variant_sets": [{"type": "DEL", "number": 1,
                                           "length_ranges": [[5, 5]],
                                           "blacklist_region_type": "ALR"}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {
                "chr21": "ACTAATCTCTTCTCTCTTCTCTCTCCGT"}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file},
                         "blacklist_regions": ["test/inputs/example_avoid_intervals_2.bed",
                                               "test/inputs/example_avoid_intervals_3.bed",
                                               "test/inputs/example_avoid_interval.vcf"],
                         "variant_sets": [{"type": "DEL", "number": 1,
                                           "length_ranges": [[5, 5]],
                                           "blacklist_region_type": "Alu"},
                                          {"type": "DEL", "number": 2,
                                           "length_ranges": [[5, 5]],
                                           "blacklist_region_type": "L1"}
                                          ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {
                "chr21": "ACTAATCTCTTCTCTCTTCTCTCTCCGT"}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file},
                         # can specify same files for blacklist and overlap regions
                         "overlap_regions": ["test/inputs/example_avoid_intervals_2.bed"],
                         "blacklist_regions": ["test/inputs/example_avoid_intervals_2.bed"],
                         "variant_sets": [{"type": "DEL", "number": 1,
                                           "length_ranges": [[1, 5]],
                                           "overlap_mode": "exact",
                                           "overlap_region_type": "L1"},
                                          {"type": "DEL", "number": 1,
                                           "length_ranges": [[5, 10]],
                                           "blacklist_region_type": "L1"}
                                          ]}],
                       self.hap1, self.hap2, self.bed)
        ]
        # ----------- test objects for unbounded dispersions -----------
        self.test_objects_unbounded_disp = [
            TestObject([self.ref_file, {"chr21": utils.generate_seq(1000)}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "variant_sets": [
                                       {"type": "TRA_NONRECIPROCAL", "number": 1,
                                        "length_ranges": [[2, 2], [1, None]]},
                                       {"type": "dDUP", "number": 1,
                                        "length_ranges": [[2, 2], [1, None]]},
                                       {"type": "INV_dDUP", "number": 1,
                                        "length_ranges": [[2, 2], [1, None]]},
                                       {"type": "dDUP_iDEL", "number": 1,
                                        "length_ranges": [[2, 2], [2, 2], [1, None]]},
                                       {"type": "INS_iDEL", "number": 1,
                                        "length_ranges": [[2, 2], [2, 2], [1, None]]}
                                   ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {
                "chr21": utils.generate_seq(1000)}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "variant_sets": [
                                       {"type": "Custom", "source": ["A", "_", "B", "_", "C"],
                                        "target": ["a", "_", "b", "b^", "_", "C"], "number": 2,
                                        "length_ranges": [[2, 2], [2, 2], [2, 2], [2, None], [2, None]]},
                                   ]}],
                       self.hap1, self.hap2, self.bed),
        ]
        # ---------- test objects with custom SVs ----------
        self.test_objects_no_dis_custom = [
            TestObject([self.ref_file, {"chr21": "AGACT"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "variant_sets": [
                                       {"type": "Custom", "source": ["A"], "target": ["a", "B"],
                                        "number": 1, "length_ranges": [[5, 5], [5, 5]]}
                                   ]}],
                       self.hap1, self.hap2, self.bed),
        ]

        self.test_objects_custom_overlap = [
            TestObject([self.ref_file, {"chr21": "GCAGACTGAC"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed,
                                   "variant_sets": [
                                       {"type": "Custom", "source": ["A", "B"], "target": ["A", "A^"], "number": 1,
                                        "length_ranges": [[5, 5], [5, 5]],
                                        "overlap_mode": "exact",
                                        "overlap_region_type": [None, "L1HS"]}
                                   ]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"chr21": "GCAGACTGAC"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": self.test_overlap_bed,
                                   "variant_sets": [
                                       {"type": "Custom", "source": ["A", "B"], "target": ["A", "A^"], "number": 1,
                                        "length_ranges": [[6, 6], [5, 5]],
                                        "overlap_mode": "exact",
                                        "overlap_region_type": [None, "L1HS"]}
                                   ]}],
                       self.hap1, self.hap2, self.bed),
        ]

        # test objects for config files specifying a minimum inter-SV breakend distance
        self.test_objects_intersv_distance = [
            TestObject([self.ref_file, {"Chr21": utils.generate_seq(600)}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file, "min_intersv_dist": 10},
                         "variant_sets": [
                             {"type": "DEL", "number": 20, "length_ranges": [(5, 5)]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"Chr21": utils.generate_seq(1000)}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file, "min_intersv_dist": 10},
                         "variant_sets": [
                             {"type": "dDUP", "number": 20, "length_ranges": [(5, 5), (5, 50)]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"Chr21": utils.generate_seq(20)}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file, "min_intersv_dist": 10},
                         "variant_sets": [
                             {"type": "DUP", "number": 3, "length_ranges": [(5, 5)]}]}],
                       self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {"Chr21": utils.generate_seq(10)}],
                       [self.par,
                        {"sim_settings": {"reference": self.ref_file},
                         "variant_sets": [
                             {"type": "DUP", "number": 2, "length_ranges": [(5, 5)]}]}],
                       self.hap1, self.hap2, self.bed)
            ]

        self.test_objects_roi_placement_failure = [
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [10, 10],
                                                     "overlap_mode": "exact"}]}],
                       self.hap1, self.hap2, self.bed),
            # combine two input files, filter all by length
            TestObject([self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed,
                                                       self.test_overlap_bed_2],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [10, 10],
                                                     "overlap_mode": "exact"}]}],
                       self.hap1, self.hap2, self.bed),
            # combine two input files, filter all by chromosome
            TestObject([self.ref_file, {"chr19": "CTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed,
                                                       self.test_overlap_bed_2],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [10, 10],
                                                     "overlap_mode": "exact"}]}],
                       self.hap1, self.hap2, self.bed),
            # type-specific num_overlap params
            TestObject(
                [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                [self.par,
                 {"sim_settings": {"reference": self.ref_file},
                  "overlap_regions": [self.test_overlap_bed, self.test_overlap_bed_2],
                  "variant_sets": [{"type": "DEL", "number": 1,
                                    "length_ranges": [[1, 10]]},
                                   {"type": "DEL", "number": 2,
                                    "overlap_region_length_range": [1, 10],
                                    "overlap_mode": "exact",
                                    "overlap_region_type": ["L1HS"]},
                                   {"type": "DEL", "number": 4,
                                    "overlap_region_length_range": [1, 10],
                                    "overlap_mode": "exact",
                                    "overlap_region_type": ["ALR/Alpha"]}
                                   ]}],
                self.hap1, self.hap2, self.bed),
            # type-specific num_overlap param > num available (ALR)
            TestObject(
                [self.ref_file, {"chr21": "CTCCGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                [self.par,
                 {"sim_settings": {"reference": self.ref_file},
                  "overlap_regions": [self.test_overlap_bed, self.test_overlap_bed_2],
                  "variant_sets": [{"type": "DEL", "number": 2,
                                    "overlap_region_length_range": [1, 5],
                                    "overlap_mode": "exact",
                                    "overlap_region_type": ["L1HS"]},
                                   {"type": "DEL", "number": 3,
                                    "overlap_region_length_range": [1, 5],
                                    "overlap_mode": "exact",
                                    "overlap_region_type": ["ALR/Alpha"]}
                                   ]}],
                self.hap1, self.hap2, self.bed),
            TestObject([self.ref_file, {
                "chr21": "CCTCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTATCCGTCGTACTAAGTCGTACTAAGTCGTACTCCGTCGTACTAAGTCGTA"}],
                       [self.par, {"sim_settings": {"reference": self.ref_file},
                                   "overlap_regions": [self.test_overlap_bed_11],
                                   "variant_sets": [{"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [2, 4],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["Alu"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [2, 4],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["L1"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [2, 4],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["L2"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [2, 4],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["SVA"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [2, 4],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["HERVK"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [6, 8],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["Alu"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [6, 8],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["L1"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [6, 8],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["L2"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [6, 8],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["SVA"]},
                                                    {"type": "DEL", "number": 1,
                                                     "overlap_region_length_range": [6, 8],
                                                     "overlap_mode": "exact",
                                                     "overlap_region_type": ["HERVK"]}
                                                    ]}],
                       self.hap1, self.hap2, self.bed)
        ]

        self.simple_test_data = [
            ["TC", {"type": "INS", "length_ranges": [[1, 1]]}, ["TTC", "TCC", "TGC", "TAC"]],

            ["TCGA", {"type": "DEL", "length_ranges": [[2, 2]]}, ["TC", "TA", "GA"]],
            [["TCGA", "CTG"], {"type": "DEL", "length_ranges": [[2, 2]]},
             [("TC", "CTG"), ("TA", "CTG"), ("GA", "CTG"), ("TCGA", "C"), ("TCGA", "G")]],
            ["TCGA", {"type": "DEL", "length_ranges": [[3, 3]]}, ["T", "A"]],

            ["TC",   {"type": "SNP", "length_ranges": [[1, 1]]}, ["TT", "TG", "TA", "AC", "GC", "CC"]],
            ["T",    {"type": "SNP", "length_ranges": [[1, 1]]}, ["C", "G", "A"]],

            ["TCGA", {"type": "INV", "length_ranges": [[3, 4]]}, ["CGAA", "TTCG", "TCGA"]],

            ["TCGA", {"type": "DUP", "length_ranges": [[3, 4]]}, ["TCGTCGA", "TCGACGA", "TCGATCGA"]],

            ["TCGA", {"type": "TRA_NONRECIPROCAL", "length_ranges": [[2, 2], [1, None]]},
             ["GTCA", "TGAC", "GATC", "CGTA", "TACG"]],

            ["TCGA", {"type": "TRA_RECIPROCAL", "length_ranges": [[2, 2], [2, 2], [None, None]]},
             ["GATC"]],
            ["TCGA", {"type": "TRA_RECIPROCAL", "length_ranges": [[2, 2], ["A", "A"], [None, None]]},
             ["GATC"]],
            ["TCGA", {"type": "TRA_RECIPROCAL", "length_ranges": [["B", "B"], [2, 2], [None, None]]},
             ["GATC"]],
            ["TCGA", {"type": "TRA_RECIPROCAL", "length_ranges": [[3, 3], ["A-2", "A-2"], [None, None]]},
             ["ATCG", "CGAT"]],
            ["TCGA", {"type": "TRA_RECIPROCAL", "length_ranges": [[1, 1], [1, 1], [2, None]]},
             ["ACGT"]],

            ["TCGA", {"type": "dupINVdup", "length_ranges": [[1, 1], [2, 2], [1, 1]]},
             ["TTCGAA"]],

            ["TCGA", {"type": "delINVdel", "length_ranges": [[1, 1], [2, 2], ["A", "A"]]},
             ["CG"]],

            ["TCGA", {"type": "delINVdup", "length_ranges": [[1, 1], [2, 2], [1, 1]]},
             ["TCGA"]],

            ["TCGA", {"type": "dupINVdel", "length_ranges": [[1, 1], [2, 2], [1, 1]]},
             ["TCGA"]],

            ["TCGA", {"type": "delINV", "length_ranges": [[2, 2], [1, 1]]},
             ["CA", "TT"]],

            ["TCGA", {"type": "INVdel", "length_ranges": [[2, 2], [1, 1]]},
             ["GAA", "TCG"]],

            ["TCGA", {"type": "INS_iDEL", "length_ranges": [[2, 2], [1, 1], [1, None]]},
             ["GTC", "GAC"]],

            ["TCGA", {"type": "INV_DUP3", "length_ranges": [[2, 2]]}, ["GAGAGA", "TCTCTC", "TCGCGA"]],
            ["TCGA", {"type": "INV_DUP3", "length_ranges": [[1, 1]]}, ["AACGA", "TGGGA", "TCCCA", "TCGTT"]],

            ["TCGA", {"type": "INV_DUP", "length_ranges": [[3, 3]]}, ["TCGCGAA", "TCGATCG"]],

            ["TCGA", {"type": "fldup_INV", "length_ranges": [[2, 2], [2, 2]]}, ["TCTCGA"]],

            ["TCGA", {"type": "INV_fldup", "length_ranges": [[2, 2], [2, 2]]}, ["TCGAGA"]],

            [["TCGA", "CAT"], {"type": "dDUP", "length_ranges": [[3, 3], [1, None]]},
             [("TCGATCG", "CAT"), ("CGATCGA", "CAT")]],
            [["TCGA", "CAT"], {"type": "dDUP", "length_ranges": [[3, 3], [None, None]],
                                "interchromosomal": True},
             [("TCGA"[:i] + "CAT" + "TCGA"[i:], "CAT") for i in range(1, len("TCGA"))] +
             [("TCGA", "CAT"[:j] + s + "CAT"[j:]) for j in range(1, len("CAT")) for s in ("TCG", "CGA")]],

            [["TCGA", "CAT"], {"type": "INV_dDUP", "length_ranges": [[3, 3], [1, None]]},
             [("TCGACGA", "CAT"), ("TCGTCGA", "CAT")]],

            ["TCGA", {"type": "INV_TRA", "length_ranges": [[3, 3], [1, 1]]}, ["ACGA", "TCGT"]],

            ["T", {"type": "Custom", "source": "A", "target": "AA*^", "length_ranges": [[1, 1]],
                   "divergence_prob": 1.0},
             ["TC", "TG", "TA"]],
            [["GGCCTT", "CA"], [{"type": "Custom", "source": "ABC", "target": "AC",
                               "length_ranges": [[2,2],[2,2],[2,2]]},],
             [("GGTT", "CA")]],

            [["GGTT", "CA"], {"type": "TRA_CHROM_BALANCED"},
             [("CGTT", "GA"), ("CTT", "GGA"), ("GGTG", "AA")]],

            [["GGTT", "CA"], {"type": "TRA_CHROM_UNBALANCED"},
             [("CGTT", "CA"), ("GGTT", "GA"),
              ("CTT", "CA"), ("GGTT", "GGA"),
              ("GGTG", "CA"), ("GGTT", "AA")
              ]],

            [["GGCCTT", "CA"], [{"type": "Custom", "source": "A_B_C", "target": "A__C",
                               "length_ranges": [[1,1],[2,2],[1,1],[1,1],[1,1]]},
                                {"type": "TRA_CHROM_BALANCED"}],
             [("CGTT", "GA"), ("CTT", "GGA"), ("GGG", "AAA"), ("GGTG", "AA")]],

            ["TCGA", {"type": "Custom", "source": "ABCD", "target": "c^bda^", "length_ranges": [[1,1],[1,1],[1,1],[1,1]]},
             ["CGTA"]],
        ]

    def tearDown(self):
        try:
            shutil.rmtree(self.test_dir)
        except Exception as exc:
            print(f'Error removing test dir {self.test_dir}: {exc}')

    # helper method for tests where the output will be in a known list of possibilities
    def helper_test_known_output_svs(self, config_event_obj, target_frags=None):
        # target_frags: optional input frags to be checked for match with output frags
        config = config_event_obj
        config.initialize_files()
        curr_sim = SVSimulator(config.par)
        curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
        changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
        config.remove_test_files()
        if target_frags is not None:
            if changed_frag_1 == changed_frag_2:
                self.assertIn(changed_frag_1, target_frags)
            else:
                if not (changed_frag_1 in target_frags or changed_frag_2 in target_frags):
                    print(f'{changed_frag_1=} {changed_frag_2=} {target_frags=}', file=sys.stderr)
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
        # TRA_NONRECIPROCAL -- ref: CT
        # same output for forward and backward TRA_NONRECIPROCAL
        self.helper_test_known_output_svs(self.test_dispersion_objects[0], ['TC'])
        # dDUP
        # different output for forward and backward dDUP
        self.helper_test_known_output_svs(self.test_dispersion_objects[1], ['CTC', 'TCT'])
        # # INV_dDUP
        self.helper_test_known_output_svs(self.test_dispersion_objects[2], ['CTG', 'ACT'])
        # # dDUP_iDEL
        self.helper_test_known_output_svs(self.test_dispersion_objects[3], ['CTC', 'GTG'])
        # # INS_iDEL
        self.helper_test_known_output_svs(self.test_dispersion_objects[4], ['TC', 'GT'])
        # # reciprocal TRA_RECIPROCAL -- ref: CTTTA
        self.helper_test_known_output_svs(self.test_dispersion_objects[5], ['ATTTC'])
        # # reciprocal TRA_RECIPROCAL, unbounded above, bounded below -- ref: CTTTA
        self.helper_test_known_output_svs(self.test_dispersion_objects[6], ['ATTTC'])

        with self.assertRaises(Exception):
            self.helper_test_known_output_svs(self.test_dispersion_objects[7], ['TATCT'])

        self.helper_test_known_output_svs(self.test_dispersion_objects[8], ['TATCT'])
        self.helper_test_known_output_svs(self.test_dispersion_objects[9], ['TACTT', 'TTACT'])
        self.helper_test_known_output_svs(self.test_dispersion_objects[10], ['TATCT', 'TTCTA', 'CTATT'])

        self.helper_test_known_output_svs(self.test_dispersion_objects[11], ['TATCT', 'TTCTA', 'CTATT'])

    def test_overlap_placement_simple(self):
        # simple events
        for i, config in enumerate(self.test_objects_overlap_simple):
            # -- testing new overlap flag on first test object --
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
            if i == 0:
                self.assertTrue('CTGTCGTA' in [changed_frag_1, changed_frag_2])
            elif i == 1:
                #self.assertEqual(len(curr_sim.overlap_regions.roi_dict.values()), 1)
                self.assertTrue('CCCC' in changed_frag_1 or 'CCCC' in changed_frag_2)
            elif i == 2:
                self.assertTrue('AA' not in changed_frag_1 or 'AA' not in changed_frag_2)


    def test_overlap_placement_complex(self):
        # complex events
        for i, config in enumerate(self.test_objects_overlap_cplx):
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)

            def sv_source_segments(sv):
                return sorted([(operation.source_region.start, operation.source_region.end)
                               for operation in sv.operations])

            def sv_target_segments(sv):
                return sorted([(operation.target_region.start, operation.target_region.end)
                               for operation in sv.operations])

            if i == 0:
                self.helper_test_known_output_svs(self.test_objects_overlap_cplx[i], ['CTGATGA', 'CGATGAT'])
            elif i == 1:
                # TRA_NONRECIPROCAL [5,10) -> [10]; INV_dDUP [1,2) -> [0] or [3] - source ref: CTGATATGGAC
                changed_frag_1, changed_frag_2 = self.helper_test_known_output_svs(self.test_objects_overlap_cplx[i])
                # need to account for the events being placed on opposite haplotypes, so will check for each separately
                # --> check for INV_dDUP in first four characters of output refs
                self.assertTrue(changed_frag_1[:4] in ['CTGA', 'ACTG'] or changed_frag_2[:4] in ['CTGA', 'ACTG'])
                # --> check for TRA_NONRECIPROCAL in second half of refs
                self.assertTrue(changed_frag_1[-7:] in ['TCATGGA', 'ATGGATC'] or changed_frag_2[-7:] in ['TCATGGA', 'ATGGATC'])
            elif i == 2:
                for sv in curr_sim.svs:
                    self.assertIn(sv_source_segments(sv),
                                  [[(10 + i * 3, 13 + i * 3), (13 + i * 3, 16 + i * 3), (16 + i * 3, 19 + i * 3)] for i
                                   in range(3)])
            elif i == 3:
                for sv in curr_sim.svs:
                    if sv.info['SVTYPE'] == 'delINV':
                        self.assertIn(sv_source_segments(sv),
                                      [[(15 + i * 3, 18 + i * 3), (18 + i * 3, 21 + i * 3)] for i in range(3)])
                    elif sv.info['SVTYPE'] == 'INVdel':
                        self.assertIn(sv_source_segments(sv),
                                      [[(8 + i * 2, 10 + i * 2), (10 + i * 2, 12 + i * 2)] for i in range(3)])
            elif i == 4:
                # partial overlap for flanked INV
                possible_shifts = [(10 + i * 3, 13 + i * 3, 16 + i * 3) for i in range(3)]
                correct_shift = False
                (frag_a, frag_b, frag_c) = tuple(seg[0] for seg in sv_source_segments(curr_sim.svs[0]))
                for (ivl_a, ivl_b, ivl_c) in possible_shifts:
                    if all(np.abs(frag - ivl) < 3 for (frag, ivl) in
                           zip((ivl_a, ivl_b, ivl_c), (frag_a, frag_b, frag_c))):
                        correct_shift = True
                self.assertTrue(correct_shift)
            elif i in [5, 6, 7]:
                dispersion_target = [
                    operation.target_region
                    for operation in curr_sim.svs[0].operations
                    if not operation.transform.is_in_place][0]
                self.assertTrue(13 <= dispersion_target.start <= 16)
            elif i == 8:
                inv_ddup = [sv for sv in curr_sim.svs if sv.info['SVTYPE'] == 'INV_dDUP'][0]
                dispersion_target = [
                    operation.target_region
                    for operation in inv_ddup.operations
                    if not operation.transform.is_in_place][0]
                self.assertTrue(10 <= dispersion_target.start <= 13 or
                                16 <= dispersion_target.start <= 20)
            elif i == 9:
                # flanked-INV with anchor = full SV
                self.assertEqual(sv_source_segments([sv for sv in curr_sim.svs
                                                     if sv.info['SVTYPE'] == 'delINVdel'][0]),
                                 [(12, 14), (14, 17), (17, 19)])
            elif i == 10:
                self.assertEqual(sv_source_segments(curr_sim.svs[0]),
                                 [(12, 15), (15, 19)])
            elif i == 11:
                frag_bounds = sv_source_segments(sv)
                self.assertTrue(is_overlapping([(12, 19)], (frag_bounds[0][0], frag_bounds[-1][1])))

    def test_frag_level_overlap(self):
        for i, config in enumerate(self.test_objects_frag_level_overlap):
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            def sym_src(sym):
                return [operation.source_region
                        for operation in curr_sim.svs[0].operations
                        if operation.op_info['symbol'].name == sym][0]
            if i <= 10:
                src_ev = sym_src('A')
                disp_ev = utils.Region(chrom='chr21', start=curr_sim.svs[0].placement[1].pos,
                                       end=curr_sim.svs[0].placement[2].pos)
                #disp_ev = [operation for operation in curr_sim.svs[0].operations
                #           if operation.op_info['symbol'].name.startswith('_')][0]

            # exact source match for dDUP, INV_dDUP, TRA_NONRECIPROCAL
            if i in [0, 1, 2]:
                self.assertEqual((2, 4), (src_ev.start, src_ev.end))
            # exact source match for dDUP, INV_dDUP, TRA_NONRECIPROCAL
            if i in [3, 4]:
                b_ev = sym_src('B')
                self.assertEqual((13, 15), ((src_ev.start, src_ev.end) if i == 3 else (b_ev.start, b_ev.end)))
            # partial source match for dDUP, INV_dDUP, TRA_NONRECIPROCAL
            if i in [5, 6, 7]:
                self.assertTrue(is_overlapping([(2, 4), (5, 10)], (src_ev.start, src_ev.end)))
            # target match for dDUP, INV_dDUP, TRA_NONRECIPROCAL
            if i in [8, 9, 10]:
                self.assertTrue(is_overlapping([(2, 5), (5, 11), (11, 13)], (disp_ev.start, disp_ev.end)))
            # delINVdel, delINV, INVdel: specify any single fragment
            if i in [11, 12, 13, 14, 15]:
                ovlp_ev = None
                if i in [11, 15]:
                    ovlp_ev = sym_src('B')
                elif i in [12, 14]:
                    ovlp_ev = sym_src('A')
                elif i == 13:
                    ovlp_ev = sym_src('C')
                self.assertEqual((5, 10), (ovlp_ev.start, ovlp_ev.end))
            if i == 16:
                b_ev = sym_src('B')
                self.assertEqual((5, 10), (b_ev.start, b_ev.end))

    def test_partial_overlap_placement(self):
        for i in range(len(self.test_objects_partial_overlap)):
            if i not in [0, 1]: continue
            config = self.test_objects_partial_overlap[i]
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
            # to simplify checking the different possible outcomes of each test: sorting SVs by start position
            curr_sim.svs.sort(key=lambda x: get_span(x).start)
            changed_frag_1, changed_frag_2 = config.get_actual_frag(return_haps='both')
            if i == 0:
                # source: CTCCGTAGTA -> four possible valid del intervals, checking each
                self.assertTrue(any(possible_del not in changed_frag_1 for possible_del in ['CT', 'CC', 'CG', 'GT']) or
                                any(possible_del not in changed_frag_2 for possible_del in ['CT', 'CC', 'CG', 'GT']))
            if i == 1:
                # check that between the two SVs one is a perfect match and is other overlaps
                case_a = (get_span(curr_sim.svs[0]).start, get_span(curr_sim.svs[0]).end) == (2, 4) and is_overlapping([(15, 20)], (
                get_span(curr_sim.svs[1]).start, get_span(curr_sim.svs[1]).end))
                case_b = (get_span(curr_sim.svs[1]).start, get_span(curr_sim.svs[1]).end) == (15, 20) and is_overlapping([(2, 4)], (
                get_span(curr_sim.svs[0]).start, get_span(curr_sim.svs[0]).end))
                self.assertTrue(case_a or case_b)
            if i == 2:
                case_a = (get_span(curr_sim.svs[0]).start, get_span(curr_sim.svs[0]).end) == (2, 4) and is_overlapping([(22, 25)], (
                get_span(curr_sim.svs[2]).start, get_span(curr_sim.svs[2]).end))
                case_b = (get_span(curr_sim.svs[2]).start, get_span(curr_sim.svs[2]).end) == (22, 25) and is_overlapping([(2, 4)], (
                get_span(curr_sim.svs[0]).start, get_span(curr_sim.svs[0]).end))
                self.assertTrue(case_a or case_b)
                self.assertTrue((get_span(curr_sim.svs[1]).start, get_span(curr_sim.svs[1]).end) == (11, 15))
            if i == 3:
                sv_intervals = [(get_span(sv).start, get_span(sv).end) for sv in curr_sim.svs]
                self.assertTrue((11, 15) in sv_intervals)
                partial_ovl_sv = curr_sim.svs[(sv_intervals.index((11, 15)) + 1) % 2]
                self.assertTrue(is_overlapping([(2, 4)], (get_span(partial_ovl_sv).start, get_span(partial_ovl_sv).end)) or
                                is_overlapping([(22, 25)], (get_span(partial_ovl_sv).start, get_span(partial_ovl_sv).end)))

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
        # test of avoiding pre-specified intervals
        config_no_dis = self.test_objects_no_dis[0]
        config_no_dis.initialize_files()
        curr_sim = SVSimulator(config_no_dis.par)
        curr_sim.produce_variant_genome(config_no_dis.hap1,
                                        config_no_dis.hap2, config_no_dis.ref, config_no_dis.bed, export_to_file=False)

        def sv_source_segments(sv):
            return sorted([(operation.source_region.start, operation.source_region.end)
                           for operation in sv.operations])

        source_segs = sv_source_segments(curr_sim.svs[0])
        self.assertTrue(all(0 <= seg_start < seg_end <= 16
                            for seg_start, seg_end in source_segs))

    def test_req_space_filtering(self):
        config_req_space = self.test_objects_req_space[0]
        config_req_space.initialize_files()
        curr_sim = SVSimulator(config_req_space.par)
        with self.assertRaises(Exception):
            curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta)

    def test_unbounded_dispersion(self):
        for i, config in enumerate(self.test_objects_unbounded_disp):
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            #curr_sim.sim_settings['max_tries'] = 2000
            curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)


    def test_roi_placement_failure(self):
        for config in self.test_objects_roi_placement_failure:
            config.initialize_files()
            curr_sim = SVSimulator(config.par)
            with self.assertRaises(Exception):
                curr_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)

    def test_simple_tests(self):
        for test_num, (ref, vs_config, expected_outputs) in enumerate(self.simple_test_data):
            if not isinstance(vs_config, list):
                vs_config = [vs_config]

            test_object = TestObject([self.ref_file, {f"chrTest{ref_num}": ref_seq
                                                      for ref_num, ref_seq in enumerate(as_list(ref))}],
                                     [self.par, {"sim_settings": {"reference": self.ref_file,
                                                                  "homozygous_only": True,
                                                                  "min_intersv_dist": 0},
                                                 "variant_sets": [dict(number=1, **vs_conf)
                                                                  for vs_conf in vs_config]}],
                                     self.hap1, self.hap2, self.bed)

            results_seen = set()
            attempt_num = 0
            expected_results = set(expected_outputs)
            while len(results_seen) < len(expected_results) and attempt_num < len(expected_results) * 100:
                attempt_num += 1
                results = self.helper_test_known_output_svs(
                    test_object, expected_results)
                results_seen.update(results)
            assert results_seen == expected_results, f'{test_num=} {vs_config=}'

def test_inv(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "a.yaml"
    p.write_text("""
sim_settings:
    reference: "test/inputs/test01.fa"
    max_tries: 100
variant_sets:
    - type: "INV"  # "A" -> ""
      number: 1
      length_ranges: [[3, 3]]
    """)

    simulator = SVSimulator(config_path=str(p))
    simulator.run()

    svs = simulator.svs
    assert len(svs) == 1
    assert len(svs[0].operations) == 1
    assert (svs[0].operations[0].transform == 
            Transform(transform_type=TransformType.INV, is_in_place=True, divergence_prob=0))
    assert (svs[0].operations[0].target_region.replace(order_key=()) == 
            svs[0].operations[0].source_region)

def test_inv_exact(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    a_bed = d / "a.bed"
    a_bed.write_text("""
    chr19\t0\t3\tLINE1
    """)
    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test01.fa"
    max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"  # "A" -> ""
      number: 1
      length_ranges: [[null, null]]
      overlap_region_type: ["LINE1"]
      overlap_mode: exact
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs
    assert len(svs) == 1
    assert len(svs[0].operations) == 1
    assert (svs[0].operations[0].transform == 
            Transform(transform_type=TransformType.INV, is_in_place=True, divergence_prob=0))
    assert (svs[0].operations[0].target_region.replace(order_key=()) == 
            svs[0].operations[0].source_region)

    assert (svs[0].operations[0].target_region.replace(order_key=()) ==
            Region(chrom='chr19', start=0, end=3))


    
def test_inv_contained(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    a_bed = d / "a.bed"
    a_bed.write_text("""
    chr19\t3\t7\tLINE1
    """)
    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test01.fa"
    max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"  # "A" -> ""
      number: 1
      length_ranges: [[3, 3]]
      overlap_region_type: ["LINE1"]
      overlap_mode: contained
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs
    assert len(svs) == 1
    assert len(svs[0].operations) == 1
    assert (svs[0].operations[0].transform == 
            Transform(transform_type=TransformType.INV, is_in_place=True, divergence_prob=0))
    assert (svs[0].operations[0].target_region.replace(order_key=()) == 
            svs[0].operations[0].source_region)

    assert (svs[0].operations[0].target_region.replace(order_key=()) in
            [Region(chrom='chr19', start=3, end=6), Region(chrom='chr19', start=4, end=7)])

    
def test_inv_partial(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    a_bed = d / "a.bed"
    a_bed.write_text("""
    chr19\t3\t7\tLINE1
    """)
    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test01.fa"
    max_tries: 1
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "INV"  # "A" -> ""
      number: 1
      length_ranges: [[3, 3]]
      overlap_region_type: ["LINE1"]
      overlap_mode: partial
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs
    assert len(svs) == 1
    assert len(svs[0].operations) == 1
    assert (svs[0].operations[0].transform == 
            Transform(transform_type=TransformType.INV, is_in_place=True, divergence_prob=0))
    assert (svs[0].operations[0].target_region.replace(order_key=()) == 
            svs[0].operations[0].source_region)

    assert (svs[0].operations[0].target_region.replace(order_key=()) in
            [Region(chrom='chr19', start=0, end=3),
             Region(chrom='chr19', start=1, end=4),
             Region(chrom='chr19', start=2, end=5),
             Region(chrom='chr19', start=5, end=8),
             Region(chrom='chr19', start=6, end=9),])

def test_trEXP(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    a_bed = d / "a.bed"
    a_bed.write_text("""
    chrA\t4\t16\tALU\t3
    """)
    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test_tr.fa"
    max_tries: 1
    homozygous_only: true
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "trEXP"
      number: 1
      repeat_count_change_range: [2, 2]
      overlap_region_type: ["ALU"]
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs

    sim_fa = d / "sim.hapA.fa"

    assert len(svs) == 1
    with FastaFile(sim_fa) as fasta_file:
        hap = fasta_file.fetch(fasta_file.references[0])
    assert hap == 'AAAATCGTCGTCGTCGTCGTCGAAAA'
    

def test_trCON(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    a_bed = d / "a.bed"
    a_bed.write_text("""
    chrA\t4\t16\tALU\t3
    """)
    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test_tr.fa"
    max_tries: 1
    homozygous_only: true
overlap_regions: ["{a_bed}"]
variant_sets:
    - type: "trCON"
      number: 1
      repeat_count_change_range: [2, 2]
      overlap_region_type: ["ALU"]
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs

    sim_fa = d / "sim.hapA.fa"

    assert len(svs) == 1
    with FastaFile(sim_fa) as fasta_file:
        hap = fasta_file.fetch(fasta_file.references[0])
    assert hap == 'AAAATCGTCGAAAA'

    
def test_trINS(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()

    cfg = d / "a.yaml"
    cfg.write_text(f"""
sim_settings:
    reference: "test/inputs/test_tr.fa"
    max_tries: 1
    homozygous_only: true
variant_sets:
    - type: "trINS"
      number: 1
      repeat_count_change_range: [4, 4]
      repeat_sequence_choices: ["CAG"]
    """)

    simulator = SVSimulator(config_path=str(cfg))
    simulator.run()

    svs = simulator.svs

    sim_fa = d / "sim.hapA.fa"

    assert len(svs) == 1
    with FastaFile(sim_fa) as fasta_file:
        hap = fasta_file.fetch(fasta_file.references[0])

    assert ('CAG' * 4) in hap
    assert ('CAG' * 5) not in hap

if __name__ == "__main__":
    unittest.main()
