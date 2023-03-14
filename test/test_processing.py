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

    def extract_bed_records(self):
        # parse bed record into dict for easy comparison
        # --> example split bed record: ['chr19', '0', '3', 'chr19', '0', '3', 'DEL', '3', '1/1', 'DEL', '1', '0']
        bed_records = []
        with open(self.bed) as f:
            for line in f:
                ln = line.split()
                bed_record = {'source_chr': ln[0], 'source_s': ln[1], 'source_e': ln[2],
                              'target_chr': ln[3], 'target_s': ln[4], 'target_e': ln[5],
                              'ev_type': ln[6], 'len': ln[7], 'zyg': ln[8], 'parent_type': ln[9],
                              'nth_sv': ln[10], 'order': ln[11]}
                bed_records.append(bed_record)
        return bed_records


class TestProcessing(unittest.TestCase):
    def setUp(self):
        # runs before every test
        self.ref_file = "test/inputs/test.fa"
        self.par = "test/inputs/par.yaml"

        self.hap1 = "test/inputs/test1.fa"
        self.hap2 = "test/inputs/test2.fa"
        self.bed = "test/inputs/out.bed"
        self.vcf = "test/inputs/out.vcf"

        self.test_objects_simple_events = {'DEL': TestProcObject([self.ref_file, {"chr19": "CTC"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "DEL", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf),
                                           'DUP': TestProcObject([self.ref_file, {"chr19": "CTC"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "DUP", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf),
                                           'INV': TestProcObject([self.ref_file, {"chr19": "CTC"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "INV", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf)}

        self.formatter = FormatterIO(self.par)

    def test_export_bedpe_simple_events(self):
        for ev_type in ['DEL', 'DUP', 'INV']:
            config = self.test_objects_simple_events[ev_type]
            config.initialize_files()
            curr_sim = SV_Simulator(config.ref, config.par)
            curr_sim.apply_transformations(FastaFile(curr_sim.ref_file))
            self.formatter.export_to_bedpe(curr_sim.svs, self.bed)
            record = config.extract_bed_records()[0]
            self.assertTrue(record['source_chr'] == record['target_chr'] == 'chr19')
            self.assertTrue(record['source_s'] == record['target_s'] == '0')
            self.assertTrue(record['source_e'] == record['target_e'] == '3')
            self.assertTrue(record['ev_type'] == record['parent_type'] == ev_type)
            self.assertTrue(record['len'] == '3')


if __name__ == "__main__":
    unittest.main()
