import unittest
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from simulate import *
from test_simulate import TestObject


class TestKnownSVs(unittest.TestCase):

    def setUp(self):
        test_vcf_simple_del = "test/inputs/example_simple_del.vcf"
        test_vcf_simple_dup = "test/inputs/example_simple_dup.vcf"
        test_vcf_simple_inv = "test/inputs/example_simple_inv.vcf"
        test_vcf_simple_ins = "test/inputs/example_simple_ins.vcf"
        test_vcf_simple_ins_no_insseq = "test/inputs/example_simple_ins_no_insseq.vcf"

        test_vcf_multidel = "test/inputs/example_multidel.vcf"
        test_vcf_del_ins = "test/inputs/example_del_ins.vcf"
        test_vcf_del_ins_del = "test/inputs/example_del_ins_del.vcf"
        test_vcf_del_dup_del = "test/inputs/example_del_dup_del.vcf"
        test_vcf_del_inv_del = "test/inputs/example_del_inv_del.vcf"
        test_vcf_dup_dup_ins = "test/inputs/example_dup_dup_ins.vcf"
        test_vcf_div_dDUP = "test/inputs/example_div_dDUP.vcf"
        test_vcf_dDUP = "test/inputs/example_dDUP.vcf"
        test_vcf_INV_dDUP = "test/inputs/example_INV_dDUP.vcf"
        test_vcf_TRA = "test/inputs/example_TRA.vcf"
        test_vcf_multi_dispersion = "test/inputs/example_multi_dispersion.vcf"

        ref_file = "test/inputs/test.fna"
        par = "test/inputs/par.yaml"
        hap1 = "test/inputs/test1.fna"
        hap2 = "test/inputs/test2.fna"
        bed = "test/inputs/out.bed"

        self.test_objects_simple_events = {'DEL': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [par, {"SVs": [{"vcf_path": test_vcf_simple_del}]}], hap1,
                                                             hap2, bed),
                                           'DUP': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [par, {"SVs": [{"vcf_path": test_vcf_simple_dup}]}], hap1,
                                                             hap2, bed),
                                           'INV': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [par, {"SVs": [{"vcf_path": test_vcf_simple_inv}]}], hap1,
                                                             hap2, bed),
                                           'INS': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                             [par, {"SVs": [{"vcf_path": test_vcf_simple_ins}]}], hap1,
                                                             hap2, bed),
                                           'INS_2': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [par, {"SVs": [{"vcf_path": test_vcf_simple_ins_no_insseq}]}],
                                                               hap1, hap2, bed),
                                           }
        self.test_objects_multievent = {'multiDEL': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                               [par, {"SVs": [{"vcf_path": test_vcf_multidel}]}], hap1,
                                                               hap2, bed),
                                        'del_ins': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                              [par, {"SVs": [{"vcf_path": test_vcf_del_ins}]}], hap1,
                                                              hap2, bed),
                                        'del_ins_del': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [par, {"SVs": [{"vcf_path": test_vcf_del_ins_del}]}],
                                                                  hap1, hap2, bed),
                                        'del_dup_del': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [par, {"SVs": [{"vcf_path": test_vcf_del_dup_del}]}],
                                                                  hap1, hap2, bed),
                                        'del_inv_del': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [par, {"SVs": [{"vcf_path": test_vcf_del_inv_del}]}],
                                                                  hap1, hap2, bed),
                                        'dup_dup_ins': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [par, {"SVs": [{"vcf_path": test_vcf_dup_dup_ins}]}],
                                                                  hap1, hap2, bed)}
        self.test_objects_dispersions = {'div_dDUP': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [par, {"SVs": [{"vcf_path": test_vcf_div_dDUP}]}],
                                                                hap1, hap2, bed),
                                         'dDUP': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                            [par, {"SVs": [{"vcf_path": test_vcf_dDUP}]}],
                                                            hap1, hap2, bed),
                                         'INV_dDUP': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                [par, {"SVs": [{"vcf_path": test_vcf_INV_dDUP}]}],
                                                                hap1, hap2, bed),
                                         'TRA': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                           [par, {"SVs": [{"vcf_path": test_vcf_TRA}]}],
                                                           hap1, hap2, bed),
                                         'multievent': TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                                  [par,
                                                                   {"SVs": [{"vcf_path": test_vcf_multi_dispersion}]}],
                                                                  hap1, hap2, bed)}

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

    def test_multiple_events(self):
        # both DELs heterozygous
        self.helper_test_simple_sv(self.test_objects_multievent['multiDEL'], ['GCTATCTCCGT'])
        # hom. DEL, het. INS
        self.helper_test_simple_sv(self.test_objects_multievent['del_ins'], ['GCGGGGGGGACTATCTCGT', 'GCACTATCTCGT'])
        # het. DELs, hom. INS
        self.helper_test_simple_sv(self.test_objects_multievent['del_ins_del'],
                                   ['GCGGGGGGGTATCTCGT', 'GCACGGGGGGGTATCTCGT',
                                    'GCGGGGGGGTATCTCTCCGT', 'GCACGGGGGGGTATCTCTCCGT'])
        self.helper_test_simple_sv(self.test_objects_multievent['del_dup_del'], ['GCTATCTATCTCGT', 'GCACTATCTATCTCGT',
                                                                                 'GCTATCTCGT', 'GCACTATCTCGT'])
        self.helper_test_simple_sv(self.test_objects_multievent['del_inv_del'], ['GCATACTCGT'])
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


if __name__ == '__main__':
    unittest.main()
