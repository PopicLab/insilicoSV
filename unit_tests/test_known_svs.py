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

    # TODO: sort out how to get for different genotypes
    #  --> the get_actual_frag() method only looks in the hap1 output...
    #  ------> I guess the answer is to just add an optional flag that'll look in both haps and we'll generalize
    #  the tests to check if both haplotypes look the way they should
    def helper_test_simple_sv(self, sv_type, config_event_obj, input_frag, target_frag):
        # template test method for simple SVs
        config_simple_del = config_event_obj
        config_simple_del.initialize_files()
        fixed_sim = SV_Simulator(config_simple_del.ref, config_simple_del.par)
        self.assertEqual(fixed_sim.produce_variant_genome(fasta1_out=config_simple_del.hap1,
                                                          fasta2_out=config_simple_del.hap2,
                                                          ins_fasta=config_simple_del.ref,
                                                          bedfile=config_simple_del.bed, verbose=True), True)
        changed_frag = config_simple_del.get_actual_frag()  # fetch fragment produced by simulator
        config_simple_del.remove_test_files()  # remove any output files and .fai files
        print(f'input frag from simple {sv_type} test = {input_frag}')
        print(f'changed frag from simple {sv_type} test = {changed_frag}')
        self.assertEqual(changed_frag, target_frag)

    def test_simple_del(self):
        self.helper_test_simple_sv('DEL', self.test_object_simple_DEL, 'GCACTATCTCTCCGT', 'GCTCCGT')

    def test_simple_dup(self):
        self.helper_test_simple_sv('DUP', self.test_object_simple_DUP, 'GCACTATCTCTCCGT', 'GCACTATCTCACTATCTCTCCGT')

    def test_simple_inv(self):
        self.helper_test_simple_sv('INV', self.test_object_simple_INV, 'GCACTATCTCTCCGT', 'GCGAGATAGTTCCGT')

    # TODO: How to specify an insertion interval in our example ins vcf?
    #  --> and how to handle INSSEQ? What's the actual form of INS events we should expect?
    # def test_simple_ins(self):
    #     self.helper_test_simple_sv('INS', self.test_object_simple_INS, 'GCACTATCTCTCCGT', 'GCAGGGGGGGCTATCTCTCCGT')

if __name__ == '__main__':
    unittest.main()
