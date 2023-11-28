import unittest
import sys
import os

from insilicosv import utils

from insilicosv.simulate import *
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
        self.test_vcf_div_dDUP = "test/inputs/example_div_dDUP.vcf"
        self.test_vcf_dDUP = "test/inputs/example_dDUP.vcf"
        self.test_vcf_INV_dDUP = "test/inputs/example_INV_dDUP.vcf"
        self.test_vcf_TRA = "test/inputs/example_TRA.vcf"
        self.test_vcf_multi_dispersion = "test/inputs/example_multi_dispersion.vcf"

        self.ref_file = "test/inputs/test.fna"
        self.par = "test/inputs/par.yaml"
        self.hap1 = "test/inputs/test1.fna"
        self.hap2 = "test/inputs/test2.fna"
        self.bed = "test/inputs/out.bed"

        self.test_objects_simple_events = {'DEL': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"SVs": [{"vcf_path": self.test_vcf_simple_del}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'DUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"SVs": [{"vcf_path": self.test_vcf_simple_dup}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INV': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"SVs": [{"vcf_path": self.test_vcf_simple_inv}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INS': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"SVs": [{"vcf_path": self.test_vcf_simple_ins}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INS_2': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [self.par, {"SVs": [{"vcf_path": self.test_vcf_simple_ins_no_insseq}]}],
                                                               self.hap1, self.hap2, self.bed),
                                           'SNP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [self.par, {"SVs": [{"vcf_path": self.test_vcf_snp}]}],
                                                             self.hap1, self.hap2, self.bed),
                                           'INVdup': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [self.par, {"SVs": [{"vcf_path": self.test_vcf_invdup}]}],
                                                                self.hap1, self.hap2, self.bed),
                                           }
        self.test_objects_multievent = {'multiDEL': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [self.par, {"SVs": [{"vcf_path": self.test_vcf_multidel}]}],
                                                               self.hap1, self.hap2, self.bed),
                                        'del_ins': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                              [self.par, {"SVs": [{"vcf_path": self.test_vcf_del_ins}]}],
                                                              self.hap1, self.hap2, self.bed),
                                        'del_ins_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"SVs": [{"vcf_path": self.test_vcf_del_ins_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'del_dup_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"SVs": [{"vcf_path": self.test_vcf_del_dup_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'del_inv_del': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"SVs": [{"vcf_path": self.test_vcf_del_inv_del}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'dup_dup_ins': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par, {"SVs": [{"vcf_path": self.test_vcf_dup_dup_ins}]}],
                                                                  self.hap1, self.hap2, self.bed),
                                        'multiDEL_multiSNP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                        [self.par, {"SVs": [{"vcf_path": self.test_vcf_multidel_multisnp}]}],
                                                                        self.hap1, self.hap2, self.bed),
                                        }
        self.test_objects_dispersions = {'div_dDUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [self.par, {"SVs": [{"vcf_path": self.test_vcf_div_dDUP}]}],
                                                                self.hap1, self.hap2, self.bed),
                                         'dDUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                            [self.par, {"SVs": [{"vcf_path": self.test_vcf_dDUP}]}],
                                                            self.hap1, self.hap2, self.bed),
                                         'INV_dDUP': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [self.par, {"SVs": [{"vcf_path": self.test_vcf_INV_dDUP}]}],
                                                                self.hap1, self.hap2, self.bed),
                                         'TRA': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                           [self.par, {"SVs": [{"vcf_path": self.test_vcf_TRA}]}],
                                                           self.hap1, self.hap2, self.bed),
                                         'multievent': TestObject([self.ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [self.par,
                                                                   {"SVs": [{"vcf_path": self.test_vcf_multi_dispersion}]}],
                                                                  self.hap1, self.hap2, self.bed)}

    def helper_test_simple_sv(self, config_event_obj, target_frags=None):
        # template test method for simple SVs
        config = config_event_obj
        config.initialize_files()
        fixed_sim = SV_Simulator(config.ref, config.par)
        fixed_sim.produce_variant_genome(config.hap1, config.hap2, config.ref, config.bed, export_to_file=False)
        changed_frag_hap1, changed_frag_hap2 = config.get_actual_frag(return_haps='both')
        config.remove_test_files()
        if target_frags is not None:
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
        self.helper_test_simple_sv(self.test_objects_simple_events['INVdup'], ['GCAGATAGTAGATAGTCTCCGT'])

    def test_multiple_events(self):
        # both DELs heterozygous
        self.helper_test_simple_sv(self.test_objects_multievent['multiDEL'], ['GCTATCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_multievent['multiDEL_multiSNP'], ['AATATCTCCGT'])
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
        self.helper_test_simple_sv(self.test_objects_dispersions['div_dDUP'], ['GCACTATCTCACTATTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['dDUP'], ['GCACTATCTCACTATTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['INV_dDUP'], ['GCACTATCTCATAGTTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['TRA'], ['GCCTCACTATTCCGT'])
        self.helper_test_simple_sv(self.test_objects_dispersions['multievent'], ['GACGGCCTACTATCACTATTCCGT'])

    def test_get_original_pos(self):
        for i in range(20):
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_simple_del, i, 'chr21') == i + (7 if i >= 2 else 0))
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_simple_dup, i, 'chr21') == i - (7 if i >= 9 else 0))
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_simple_inv, i, 'chr21') == (i if i < 2 or i >= 9 else 10 - i))
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_del_dup, i, 'chr21') == i + (3 if 1 < i < 11 else 0))
        self.assertTrue(utils.get_original_base_pos(self.test_vcf_multidel, 1, 'chr21'), 1)
        self.assertTrue(utils.get_original_base_pos(self.test_vcf_multidel, 9, 'chr21'), 11)
        self.assertTrue(utils.get_original_base_pos(self.test_vcf_multidel, 13, 'chr21'), 17)
        for i in range(20):
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_del_dup_2, i, 'chr21') == i + 2 * (i > 0) - 4 * (i > 5))
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_del_dup_del, i, 'chr21') == i + 2 * (i > 0) - 4 * (i > 5) + 3 * (i > 10))
            self.assertTrue(utils.get_original_base_pos(self.test_vcf_del_inv_del, i, 'chr21') == ((9 - i) if 1 < i < 6 else (i + 2 * (i > 0) + 3 * (i > 6))))


if __name__ == '__main__':
    unittest.main()
