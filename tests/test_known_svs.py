import os
import shutil
import sys
import tempfile
import unittest

from insilicosv import utils

from insilicosv.simulate import SVSimulator
from test_simulate import TestObject


class TestKnownSVs(unittest.TestCase):

    def setUp(self):
        self.test_vcf_simple_del = "test/inputs/example_simple_del.vcf"
        self.test_vcf_simple_dup = "test/inputs/example_simple_dup.vcf"
        self.test_vcf_simple_inv = "test/inputs/example_simple_inv.vcf"
        self.test_vcf_simple_ins = "test/inputs/example_simple_ins.vcf"
        self.test_vcf_simple_ins_no_insseq = "test/inputs/example_simple_ins_no_insseq.vcf"
        self.test_vcf_snp = "test/inputs/example_SNP.vcf"
        self.test_vcf_invdup = "test/inputs/example_invdup.vcf"

        self.test_vcf_multidel = "test/inputs/example_multidel.vcf"
        self.test_vcf_del_ins = "test/inputs/example_del_ins.vcf"
        self.test_vcf_del_ins_del = "test/inputs/example_del_ins_del.vcf"
        self.test_vcf_del_dup_del = "test/inputs/example_del_dup_del.vcf"
        self.test_vcf_del_dup = "test/inputs/example_del_dup.vcf"
        self.test_vcf_del_dup_2 = "test/inputs/example_del_dup_2.vcf"
        self.test_vcf_del_inv_del = "test/inputs/example_del_inv_del.vcf"
        self.test_vcf_dup_dup_ins = "test/inputs/example_dup_dup_ins.vcf"
        self.test_vcf_multidel_multisnp = "test/inputs/example_multiDEL_multiSNP.vcf"
        self.test_vcf_dDUP = "test/inputs/example_dDUP.vcf"
        self.test_vcf_INV_dDUP = "test/inputs/example_INV_dDUP.vcf"
        self.test_vcf_TRA_NONRECIP = "test/inputs/example_TRA_NONRECIP.vcf"
        self.test_vcf_multi_dispersion = "test/inputs/example_multi_dispersion.vcf"

        self.test_vcf_allele_freq_1 = "test/inputs/example_allele_freq_1.vcf"

        self.test_dir = tempfile.mkdtemp()
        self.ref_file = f"{self.test_dir}/test.fna"
        self.par = f"{self.test_dir}/par.yaml"

        self.hap1 = f"{self.test_dir}/test1.fna"
        self.hap2 = f"{self.test_dir}/test2.fna"
        self.bed = f"{self.test_dir}/out.bed"

        self.test_objects_simple_events = {'DEL': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"reference": self.ref_file,
                                                                         "variant_sets": [{"import": self.test_vcf_simple_del}]}], self.hap1, self.hap2, self.bed),
                                           'DUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"reference": self.ref_file,
                                                                         "variant_sets": [{"import": self.test_vcf_simple_dup}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INV': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"reference": self.ref_file,
                                                                         "variant_sets": [{"import": self.test_vcf_simple_inv}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INS': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"reference": self.ref_file,
                                                                         "variant_sets": [{"import": self.test_vcf_simple_ins}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INS_2': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [self.par, {"reference": self.ref_file,
                                                                           "variant_sets": [{"import": self.test_vcf_simple_ins_no_insseq}]}],
                                                               self.hap1, self.hap2, self.bed),
                                           'SNP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"reference": self.ref_file,
                                                                         "variant_sets": [{"import": self.test_vcf_snp}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'DUP_INV': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [self.par, {"reference": self.ref_file,
                                                                            "variant_sets": [{"import": self.test_vcf_invdup}]}],
                                                                self.hap1, self.hap2, self.bed),
                                           }
        self.test_objects_multievent = {'multiDEL': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [self.par, {"reference": self.ref_file,
                                                                           "variant_sets": [{"import": self.test_vcf_multidel}]}],
                                                               self.hap1, self.hap2, self.bed),
                                        'del_ins': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                              [self.par, {"reference": self.ref_file,
                                                                          "variant_sets": [{"import": self.test_vcf_del_ins}]}],
                                                              self.hap1, self.hap2, self.bed),
                                        'del_ins_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"reference": self.ref_file,
                                                                              "variant_sets": [{"import": self.test_vcf_del_ins_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'del_dup_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"reference": self.ref_file,
                                                                              "variant_sets": [{"import": self.test_vcf_del_dup_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'del_inv_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"reference": self.ref_file,
                                                                              "variant_sets": [{"import": self.test_vcf_del_inv_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'dup_dup_ins': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"reference": self.ref_file,
                                                                              "variant_sets": [{"import": self.test_vcf_dup_dup_ins}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'multiDEL_multiSNP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                        [self.par, {"reference": self.ref_file,
                                                                                    "variant_sets": [{"import": self.test_vcf_multidel_multisnp}]}],
                                                                        self.hap1, self.hap2, self.bed),
                                        }
        self.test_objects_dispersions = {'dDUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                            [self.par, {"reference": self.ref_file,
                                                                        "variant_sets": [{"import": self.test_vcf_dDUP}]}],
                                                            self.hap1, self.hap2, self.bed),
                                         'INV_dDUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [self.par, {"reference": self.ref_file,
                                                                            "variant_sets": [{"import": self.test_vcf_INV_dDUP}]}],
                                                                self.hap1, self.hap2, self.bed),
                                         'TRA': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                           [self.par, {"reference": self.ref_file,
                                                                       "variant_sets": [{"import": self.test_vcf_TRA_NONRECIP}]}],
                                                           self.hap1, self.hap2, self.bed),
                                         'multievent': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par,
                                                                   {"reference": self.ref_file,
                                                                    "variant_sets": [{"import": self.test_vcf_multi_dispersion}]}],
                                                                  self.hap1, self.hap2, self.bed)}

        self.test_objects_allele_freq = {'deterministic': TestObject([self.ref_file, {"chr21": "GCACTATCTCTGCACTATCTCGCACTATCTCTT"}],
                                                                     [self.par, {"reference": self.ref_file,
                                                                                 "variant_sets": [{"import": self.test_vcf_allele_freq_1}]}],
                                                                     self.hap1, self.hap2, self.bed)}

    def tearDown(self):
        try:
            shutil.rmtree(self.test_dir)
        except Exception as exc:
            print(f'Error removing test dir {self.test_dir}: {exc}')

    def helper_test_simple_sv(self, config_event_obj, target_frags=None):
        # template test method for simple SVs
        config = config_event_obj
        config.initialize_files()
        fixed_sim = SVSimulator(config.par)
        fixed_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed)
        changed_frag_hap1, changed_frag_hap2 = config.get_actual_frag(return_haps='both')
        config.remove_test_files()
        if target_frags is not None:
            print('HELPER', changed_frag_hap1, changed_frag_hap2, target_frags)
            self.assertTrue(changed_frag_hap1 in target_frags or changed_frag_hap2 in target_frags)
        return changed_frag_hap1, changed_frag_hap2

    def test_simple_events(self):
        self.helper_test_simple_sv(self.test_objects_simple_events['DEL'], ['GCCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_simple_events['DUP'], ['GCACTATCTACTATCTCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_simple_events['INV'], ['GCAGATAGTCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_simple_events['INS'], ['GCGGGGGGGACTATCTCTCCGT'])
        # insertion without specified INSSEQ (insertion sequence randomly generated)
        frag1, frag2 = self.helper_test_simple_sv(self.test_objects_simple_events['INS_2'])
        self.assertTrue(len(frag1) == 22 and len(frag2) == 22)
        self.helper_test_simple_sv(self.test_objects_simple_events['SNP'], ['GCGCTATCTCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_simple_events['DUP_INV'], ['GCAGATAGTAGATAGTCTCCGT'])

    def test_multiple_events(self):
        # both DELs heterozygous
        self.helper_test_simple_sv(self.test_objects_multievent['multiDEL'], ['GCTATCTCCGT'])
        #self.helper_test_simple_sv(self.test_objects_multievent['multiDEL_multiSNP'], ['AATATCTCCGT'])
        # hom. DEL, het. INS
        self.helper_test_simple_sv(self.test_objects_multievent['del_ins'], ['GCGGGGGGGACTATCTCGT', 'GCACTATCTCGT'])
        # het. DELs, hom. INS
        self.helper_test_simple_sv(self.test_objects_multievent['del_ins_del'],
                                   ['GCGGGGGGGTATCTCGT', 'GCACGGGGGGGTATCTCGT',
                                    'GCGGGGGGGTATCTCTCCGT', 'GCACGGGGGGGTATCTCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_multievent['del_dup_del'], ['GCTATCTATCTCGT', 'GCACTATCTATCTCGT',
                                                                                 'GCTATCTCGT', 'GCACTATCTCGT'])
        self.helper_test_simple_sv(self.test_objects_multievent['del_inv_del'], ['GCGATATCGT'])
        # het., hom., het.
        self.helper_test_simple_sv(self.test_objects_multievent['dup_dup_ins'], ['GCACACGGGGGGGTATCTCTCCTCCGT',
                                                                                 'GCACACTATCTCTCCTCCGT',
                                                                                 'GCACGGGGGTATCTCTCCTCCGT'])

    def test_dispersion_events(self):
        self.helper_test_simple_sv(self.test_objects_dispersions['dDUP'], ['GCACTATCTCACTATTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['INV_dDUP'], ['GCACTATCTCATAGTTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['TRA'], ['GCCTCACTATTCCGT'])
    #     self.helper_test_simple_sv(self.test_objects_dispersions['multievent'], ['GACGGCCTACTATCACTATTCCGT'])

    def test_allele_freq(self):
        self.helper_test_simple_sv(self.test_objects_allele_freq['deterministic'], 'GCACTATCTCTGCAGCACTATCTCTT')


if __name__ == '__main__':
    unittest.main()
