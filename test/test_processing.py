from simulate import SV_Simulator
from processing import FormatterIO
from test_simulate import TestObject
from pysam import FastaFile
import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestProcObject(TestObject):
    def __init__(self, ref, par, hap1, hap2, bed, vcf):
        self.vcf = vcf
        super().__init__(ref, par, hap1, hap2, bed)

    def extract_last_record(self, n, file='bed'):
        # helper method: extract the last n record from the output bed or vcf
        pass


class TestProcessing(unittest.TestCase):
    def setUp(self):
        # runs before every test
        self.ref_file = "test/inputs/test.fa"
        self.par = "test/inputs/par.yaml"

        self.hap1 = "test/inputs/test1.fa"
        self.hap2 = "test/inputs/test2.fa"
        self.bed = "test/inputs/out.bed"

        self.test_objects = [TestObject([self.ref_file, {"chr19": "CTC"}],
                                        [self.par, {"sim_settings": {"prioritize_top": True},
                                                    "SVs": [{"type": "DEL", "number": 1, "max_length": 3, "min_length": 3}]}],
                                        self.hap1, self.hap2, self.bed)]

        self.formatter = FormatterIO(self.par)

    def test_export_bedpe_simple_events(self):
        config = self.test_objects[0]
        config.initialize_files()
        curr_sim = SV_Simulator(config.ref, config.par)
        curr_sim.apply_transformations(FastaFile(curr_sim.ref_file))


if __name__ == "__main__":
    unittest.main()
