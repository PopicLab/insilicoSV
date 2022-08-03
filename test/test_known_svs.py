import unittest
import sys
import os

# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from simulate import *
from pysam import VariantFile
from unit_tests.test_simulate import TestObject


# -- Class attributes from TestObject --
# self.ref = ref[0]
# self.ref_content = ref[1]
# self.par = par[0]
# self.par_content = par[1]
# self.hap1 = hap1
# self.hap2 = hap2
# self.bed = bed
class TestKnownSVs(unittest.TestCase):

    def setUp(self):
        # runs before every test
        test_vcf_simple_del = "unit_tests/inputs/example_simple_del.vcf"
        test_vcf_simple_dup = "unit_tests/inputs/example_simple_dup.vcf"
        test_vcf_simple_inv = "unit_tests/inputs/example_simple_inv.vcf"
        test_vcf_simple_ins = "unit_tests/inputs/example_simple_ins.vcf"

        test_vcf_multidel = "unit_tests/inputs/example_multidel.vcf"
        test_vcf_del_ins = "unit_tests/inputs/example_del_ins.vcf"
        test_vcf_del_ins_del = "unit_tests/inputs/example_del_ins_del.vcf"
        test_vcf_del_dup_del = "unit_tests/inputs/example_del_dup_del.vcf"
        test_vcf_del_inv_del = "unit_tests/inputs/example_del_inv_del.vcf"
        test_vcf_dup_dup_ins = "unit_tests/inputs/example_dup_dup_ins.vcf"

        ref_file = "unit_tests/inputs/test.fna"
        par = "unit_tests/inputs/par.yaml"
        hap1 = "unit_tests/inputs/test1.fna"
        hap2 = "unit_tests/inputs/test2.fna"
        bed = "unit_tests/inputs/out.bed"

        self.test_object_simple_DEL = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                 [par, {"SVs": [{"vcf_path": test_vcf_simple_del}]}], hap1, hap2, bed)
        self.test_object_simple_DUP = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                 [par, {"SVs": [{"vcf_path": test_vcf_simple_dup}]}], hap1, hap2, bed)
        self.test_object_simple_INV = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                 [par, {"SVs": [{"vcf_path": test_vcf_simple_inv}]}], hap1, hap2, bed)
        self.test_object_simple_INS = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                 [par, {"SVs": [{"vcf_path": test_vcf_simple_ins}]}], hap1, hap2, bed)
        self.test_object_multiDEL = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                 [par, {"SVs": [{"vcf_path": test_vcf_multidel}]}], hap1, hap2, bed)
        self.test_object_del_ins = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                               [par, {"SVs": [{"vcf_path": test_vcf_del_ins}]}], hap1, hap2, bed)
        self.test_object_del_ins_del = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                              [par, {"SVs": [{"vcf_path": test_vcf_del_ins_del}]}], hap1, hap2, bed)
        self.test_object_del_dup_del = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                  [par, {"SVs": [{"vcf_path": test_vcf_del_dup_del}]}], hap1, hap2, bed)
        self.test_object_del_inv_del = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                  [par, {"SVs": [{"vcf_path": test_vcf_del_inv_del}]}], hap1, hap2, bed)
        self.test_object_dup_dup_ins = TestObject([ref_file, {"chr21": "GCACTATCTCTCCGT"}],
                                                  [par, {"SVs": [{"vcf_path": test_vcf_dup_dup_ins}]}], hap1, hap2, bed)

    # TODO: add tests for different genotypes
    def helper_test_simple_sv(self, sv_type, config_event_obj, input_frag, target_frag):
        # template test method for simple SVs
        config = config_event_obj
        config.initialize_files()
        fixed_sim = SV_Simulator(config.ref, config.par)
        self.assertEqual(fixed_sim.produce_variant_genome(fasta1_out=config.hap1, fasta2_out=config.hap2,
                                                          ins_fasta=config.ref, bedfile=config.bed, verbose=False), True)
        changed_frag = config.get_actual_frag()  # fetch fragment produced by simulator
        config.remove_test_files()  # remove any output files and .fai files
        # print(f'input frag from simple {sv_type} test = {input_frag}')
        # print(f'changed frag from simple {sv_type} test = {changed_frag}')
        self.assertEqual(changed_frag, target_frag)

    def test_simple_del(self):
        self.helper_test_simple_sv('DEL', self.test_object_simple_DEL, 'GCACTATCTCTCCGT', 'GCCTCCGT')

    def test_simple_dup(self):
        self.helper_test_simple_sv('DUP', self.test_object_simple_DUP, 'GCACTATCTCTCCGT', 'GCACTATCTACTATCTCTCCGT')

    def test_simple_inv(self):
        self.helper_test_simple_sv('INV', self.test_object_simple_INV, 'GCACTATCTCTCCGT', 'GCAGATAGTCTCCGT')

    def test_simple_ins(self):
        self.helper_test_simple_sv('INS', self.test_object_simple_INS, 'GCACTATCTCTCCGT', 'GCGGGGGGGACTATCTCTCCGT')

    def test_multi_del(self):
        self.helper_test_simple_sv('DEL, DEL', self.test_object_multiDEL, 'GCACTATCTCTCCGT', 'GCTATCTCCGT')

    def test_del_ins(self):
        self.helper_test_simple_sv('INS, DEL', self.test_object_del_ins, 'GCACTATCTCTCCGT', 'GCGGGGGGGACTATCTCGT')

    def test_del_ins_del(self):
        self.helper_test_simple_sv('DEL, INS, DEL', self.test_object_del_ins_del, 'GCACTATCTCTCCGT',
                                   'GCGGGGGGGTATCTCGT')

    def test_del_dup_del(self):
        self.helper_test_simple_sv('DEL, DUP, DEL', self.test_object_del_dup_del, 'GCACTATCTCTCCGT', 'GCTATCTATCTCGT')

    def test_del_inv_del(self):
        self.helper_test_simple_sv('DEL, INV, DEL', self.test_object_del_inv_del, 'GCACTATCTCTCCGT', 'GCATACTCGT')

    def test_dup_dup_inv(self):
        self.helper_test_simple_sv('DUP, DUP, INS', self.test_object_dup_dup_ins, 'GCACTATCTCTCCGT',
                                   'GCACACGGGGGGGTATCTCTCCTCCGT')

if __name__ == '__main__':
    unittest.main()
