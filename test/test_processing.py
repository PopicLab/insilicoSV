from simulate import SV_Simulator
from processing import FormatterIO
from test_simulate import TestObject
from pysam import FastaFile
from pysam import VariantFile
import unittest
import sys
import os
import utils
import constants

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
        self.ins_fasta = "test/inputs/ins_fasta.fa"

        self.test_objects_simple_events = {'DEL': TestProcObject([self.ref_file, {"chr19": "CTG"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "DEL", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf),
                                           'DUP': TestProcObject([self.ref_file, {"chr19": "CTG"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "DUP", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf),
                                           'INV': TestProcObject([self.ref_file, {"chr19": "CTG"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "INV", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf),
                                           'INS': TestProcObject([self.ref_file, {"chr19": "C"}],
                                                                 [self.par, {"sim_settings": {"prioritize_top": True},
                                                                             "SVs": [{"type": "INS", "number": 1,
                                                                                      "max_length": 3,
                                                                                      "min_length": 3}]}],
                                                                 self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_flanked_inversions = {'dupINVdup': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                            [self.par,
                                                                             {"sim_settings": {"prioritize_top": True},
                                                                              "SVs": [{"type": "dupINVdup", "number": 1,
                                                                                       "max_length": [2, 2, 2],
                                                                                       "min_length": [2, 2, 2]}]}],
                                                                            self.hap1, self.hap2, self.bed, self.vcf),
                                                'delINVdel': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                            [self.par,
                                                                             {"sim_settings": {"prioritize_top": True},
                                                                              "SVs": [{"type": "delINVdel", "number": 1,
                                                                                       "max_length": [2, 2, 2],
                                                                                       "min_length": [2, 2, 2]}]}],
                                                                            self.hap1, self.hap2, self.bed, self.vcf),
                                                'dupINVdel': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                            [self.par,
                                                                             {"sim_settings": {"prioritize_top": True},
                                                                              "SVs": [{"type": "dupINVdel", "number": 1,
                                                                                       "max_length": [2, 2, 2],
                                                                                       "min_length": [2, 2, 2]}]}],
                                                                            self.hap1, self.hap2, self.bed, self.vcf),
                                                'delINVdup': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                            [self.par,
                                                                             {"sim_settings": {"prioritize_top": True},
                                                                              "SVs": [{"type": "delINVdup", "number": 1,
                                                                                       "max_length": [2, 2, 2],
                                                                                       "min_length": [2, 2, 2]}]}],
                                                                            self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_dispersions = {'dDUP': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                [self.par,
                                                                 {"sim_settings": {"prioritize_top": True},
                                                                  "SVs": [{"type": "dDUP", "number": 1,
                                                                           "max_length": [3, 3],
                                                                           "min_length": [3, 3]}]}],
                                                                self.hap1, self.hap2, self.bed, self.vcf),
                                         'INV_dDUP': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                                    [self.par,
                                                                     {"sim_settings": {"prioritize_top": True},
                                                                      "SVs": [{"type": "INV_dDUP", "number": 1,
                                                                               "max_length": [3, 3],
                                                                               "min_length": [3, 3]}]}],
                                                                    self.hap1, self.hap2, self.bed, self.vcf),
                                         'TRA': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                               [self.par,
                                                                {"sim_settings": {"prioritize_top": True},
                                                                 "SVs": [{"type": "TRA", "number": 1,
                                                                          "max_length": [3, 3],
                                                                          "min_length": [3, 3]}]}],
                                                               self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_del_inv = {'delINV': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                              [self.par,
                                                               {"sim_settings": {"prioritize_top": True},
                                                                "SVs": [{"type": "delINV", "number": 1,
                                                                         "max_length": [3, 3],
                                                                         "min_length": [3, 3]}]}],
                                                              self.hap1, self.hap2, self.bed, self.vcf),
                                     'INVdel': TestProcObject([self.ref_file, {"chr19": "ACTGTC"}],
                                                              [self.par,
                                                               {"sim_settings": {"prioritize_top": True},
                                                                "SVs": [{"type": "INVdel", "number": 1,
                                                                         "max_length": [3, 3],
                                                                         "min_length": [3, 3]}]}],
                                                              self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_idel = {'dDUP_iDEL': TestProcObject([self.ref_file, {"chr19": "ACTGTCAG"}],
                                                              [self.par,
                                                               {"sim_settings": {"prioritize_top": True},
                                                                "SVs": [{"type": "dDUP_iDEL", "number": 1,
                                                                         "max_length": [3, 3, 2],
                                                                         "min_length": [3, 3, 2]}]}],
                                                              self.hap1, self.hap2, self.bed, self.vcf),
                                  'INS_iDEL': TestProcObject([self.ref_file, {"chr19": "ACTGTCAG"}],
                                                             [self.par,
                                                              {"sim_settings": {"prioritize_top": True},
                                                               "SVs": [{"type": "INS_iDEL", "number": 1,
                                                                        "max_length": [3, 3, 2],
                                                                        "min_length": [3, 3, 2]}]}],
                                                             self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_dup_inv = {'dup_INV': TestProcObject([self.ref_file, {"chr19": "ACTGTCAG"}],
                                                               [self.par,
                                                                {"sim_settings": {"prioritize_top": True},
                                                                 "SVs": [{"type": "dup_INV", "number": 1,
                                                                          "max_length": [4, 4],
                                                                          "min_length": [4, 4]}]}],
                                                               self.hap1, self.hap2, self.bed, self.vcf),
                                     'INV_dup': TestProcObject([self.ref_file, {"chr19": "ACTGTCAG"}],
                                                               [self.par,
                                                                {"sim_settings": {"prioritize_top": True},
                                                                 "SVs": [{"type": "INV_dup", "number": 1,
                                                                          "max_length": [4, 4],
                                                                          "min_length": [4, 4]}]}],
                                                               self.hap1, self.hap2, self.bed, self.vcf)}
        self.test_objects_INVdup = {'INVdup': TestProcObject([self.ref_file, {"chr19": "ACTG"}],
                                                             [self.par,
                                                              {"sim_settings": {"prioritize_top": True},
                                                               "SVs": [{"type": "INVdup", "number": 1,
                                                                        "max_length": 4,
                                                                        "min_length": 4}]}],
                                                             self.hap1, self.hap2, self.bed, self.vcf)}

        self.formatter = FormatterIO(self.par)

    def initialize_test(self, test_objects_dict, sv_type, output_type='bed', ins_fasta=None):
        # function to execute the shared logic for simulating SVs from test objects and generating bed/vcf output
        config = test_objects_dict[sv_type]
        config.initialize_files()
        curr_sim = SV_Simulator(config.ref, config.par)
        curr_sim.apply_transformations(FastaFile(curr_sim.ref_file))
        records = None
        if output_type == 'bed':
            if ins_fasta:
                utils.remove_file(ins_fasta)
            self.formatter.export_to_bedpe(curr_sim.svs, self.bed, ins_fasta=ins_fasta)
            records = config.extract_bed_records()
        elif output_type == 'vcf':
            pass
            # self.formatter.export_to_vcf(curr_sim, curr_sim.stats, vcffile=self.vcf)
            # TODO: write an extract_vcf_records() function
        else:
            raise ValueError('output_type must be \'bed\' or \'vcf\'')
        return records

    def extract_ins_fasta(self):
        # extracts single INSSEQ entry from self.ins_fasta (meant for use in isolated INS tests)
        # NB: chrom header output in the form ">chr{number}_{insertion_index}"
        chrom_header, insseq = None, None
        with open(self.ins_fasta) as f:
            for line in f:
                if line[0] == '>':
                    chrom_header = line
                else:
                    insseq = line.rstrip()  # <- strip trailing '\n'
        return chrom_header, insseq

    def singleton_event_bed_tests(self, records, sv_type, chrom, len):
        # helper function to perform the bed file checks shared by all singleton event tests
        self.assertTrue(all([record['source_chr'] == record['target_chr'] == chrom for record in records]))
        self.assertTrue(all([record['parent_type'] == sv_type for record in records]))
        self.assertTrue(all([record['len'] == len for record in records]))
        self.assertTrue(all([record['nth_sv'] == '1' for record in records]))
        self.assertTrue(
            all([record['order'] == str(int(record['ev_type'] in constants.NONZERO_ORDER_OPERATIONS)) for record in
                 records]))
        self.assertTrue(set([record['zyg'] for record in records]) in [{'1/1'}, {'0/1'}])

    def test_export_bedpe_simple_events(self):
        for sv_type in ['DEL', 'DUP', 'INV', 'INS']:
            record = self.initialize_test(self.test_objects_simple_events, sv_type,
                                          ins_fasta=(self.ins_fasta if sv_type == 'INS' else None))[0]
            self.assertTrue(record['ev_type'] == sv_type)
            self.singleton_event_bed_tests(records=[record], sv_type=sv_type, chrom='chr19', len='3')
            if sv_type == 'INS':
                self.assertTrue((record['source_s'], record['source_e']) in [('0', '0'), ('1', '1')])
                self.assertTrue((record['source_s'], record['source_e']) == (record['target_s'], record['target_e']))
                # test of writing of novel insertion sequences to separate output fasta
                chrom_header, insseq = self.extract_ins_fasta()
                self.assertTrue(chrom_header[1:].split('_')[0] == record['source_chr'])
                self.assertTrue(str(len(insseq)) == record['len'])
            else:
                self.assertTrue(record['source_s'] == record['target_s'] == '0')
                self.assertTrue(record['source_e'] == record['target_e'] == '3')

    def test_export_bedpe_flanked_inversions(self):
        for sv_type in ['dupINVdup', 'delINVdel', 'dupINVdel', 'delINVdup']:
            records = self.initialize_test(self.test_objects_flanked_inversions, sv_type)
            # checks of record fields that will be consistent across all four types
            self.singleton_event_bed_tests(records, sv_type, 'chr19', '2')
            for i in range(2):
                # example events set such that the source events must be placed at positions 0,2,4
                self.assertTrue(records[i]['source_s'] == str(2 * i) and records[i]['source_e'] == str(2 * i + 2))
            # type-wise checks for breakdown under fixed example event
            if sv_type == 'dupINVdup':
                self.assertTrue(records[0]['target_s'] == records[0]['target_e'] == '4')
                self.assertTrue(records[0]['ev_type'] == 'INVDUP')
                self.assertTrue(records[1]['target_s'] == '2')
                self.assertTrue(records[1]['target_e'] == '4')
                self.assertTrue(records[1]['ev_type'] == 'INV')
                self.assertTrue(records[2]['target_s'] == records[2]['target_e'] == '2')
                self.assertTrue(records[2]['ev_type'] == 'INVDUP')
            elif sv_type == 'delINVdel':
                self.assertTrue(all([record['source_s'] == record['target_s'] for record in records]))
                self.assertTrue(all([record['source_e'] == record['target_e'] for record in records]))
                self.assertTrue(records[0]['ev_type'] == 'DEL')
                self.assertTrue(records[1]['ev_type'] == 'INV')
                self.assertTrue(records[2]['ev_type'] == 'DEL')
            elif sv_type == 'dupINVdel':
                self.assertTrue(records[0]['target_s'] == records[0]['target_e'] == '4')
                self.assertTrue(all([record['source_s'] == record['target_s'] for record in records[1:]]))
                self.assertTrue(all([record['source_e'] == record['target_e'] for record in records[1:]]))
                self.assertTrue(records[0]['ev_type'] == 'INVDUP')
                self.assertTrue(records[1]['ev_type'] == 'INV')
                self.assertTrue(records[2]['ev_type'] == 'DEL')
            else:  # <- delINVdup
                self.assertTrue(all([record['source_s'] == record['target_s'] for record in records[:2]]))
                self.assertTrue(all([record['source_e'] == record['target_e'] for record in records[:2]]))
                self.assertTrue(records[2]['target_s'] == records[2]['target_e'] == '2')
                self.assertTrue(records[0]['ev_type'] == 'DEL')
                self.assertTrue(records[1]['ev_type'] == 'INV')
                self.assertTrue(records[2]['ev_type'] == 'INVDUP')

    def test_export_bedpe_dispersions(self):
        for sv_type in ['dDUP', 'INV_dDUP', 'TRA']:
            record = self.initialize_test(self.test_objects_dispersions, sv_type)[0]
            self.singleton_event_bed_tests([record], sv_type, 'chr19', '3')
            # interval checks accounting for forward or backward orientation of dispersion
            self.assertTrue((record['source_s'], record['source_e']) in [('0', '3'), ('3', '6')])
            self.assertTrue((record['target_s'], record['target_e']) in [('0', '0'), ('6', '6')])
            if sv_type == 'dDUP':
                self.assertTrue(record['ev_type'] == 'DUP')
            elif sv_type == 'INV_dDUP':
                self.assertTrue(record['ev_type'] == 'INVDUP')
            else:  # <- TRA
                self.assertTrue(record['ev_type'] == 'TRA')

    def test_export_bedpe_del_inv(self):
        for sv_type in ['delINV', 'INVdel']:
            records = self.initialize_test(self.test_objects_del_inv, sv_type)
            self.singleton_event_bed_tests(records, sv_type, 'chr19', '3')
            self.assertTrue((records[0]['source_s'], records[0]['source_e']) ==
                            (records[0]['target_s'], records[0]['target_e']) == ('0', '3'))
            self.assertTrue((records[1]['source_s'], records[1]['source_e']) ==
                            (records[1]['target_s'], records[1]['target_e']) == ('3', '6'))
            if sv_type == 'delINV':
                self.assertTrue(records[0]['ev_type'] == 'DEL')
                self.assertTrue(records[1]['ev_type'] == 'INV')
            else:  # <- INVdel
                self.assertTrue(records[0]['ev_type'] == 'INV')
                self.assertTrue(records[1]['ev_type'] == 'DEL')

    def test_export_bedpe_idel(self):
        for sv_type in ['dDUP_iDEL', 'INS_iDEL']:
            records = self.initialize_test(self.test_objects_idel, sv_type)
            self.singleton_event_bed_tests(records, sv_type, 'chr19', '3')
            # necessary locations of example event / accounting for forward or backward orientation
            self.assertTrue((records[0]['source_s'], records[0]['source_e']) == ('0', '3'))
            self.assertTrue((records[1]['source_s'], records[1]['source_e']) == ('5', '8'))
            self.assertTrue((records[0]['target_s'], records[0]['target_e']) in [('0', '3'), ('5', '5')])
            self.assertTrue((records[1]['target_s'], records[1]['target_e']) in [('3', '3'), ('5', '8')])
            for record in records:
                # the record with source == target should be the deletion
                if (record['source_s'], record['source_e']) == (record['target_s'], record['target_e']):
                    self.assertTrue(record['ev_type'] == 'DEL')
                else:
                    # for dDUP_iDEL, the source \neq target record is the DUP, for INS_iDEL it is TRA
                    if sv_type == 'dDUP_iDEL':
                        self.assertTrue(record['ev_type'] == 'DUP')
                    else:
                        self.assertTrue(record['ev_type'] == 'TRA')

    def test_export_bedpe_dup_inv(self):
        for sv_type in ['dup_INV', 'INV_dup']:
            records = self.initialize_test(self.test_objects_dup_inv, sv_type)
            self.singleton_event_bed_tests(records, sv_type, 'chr19', '4')
            self.assertTrue((records[0]['source_s'], records[0]['source_e']) == ('0', '4'))
            self.assertTrue((records[1]['source_s'], records[1]['source_e']) == ('4', '8'))
            # for dup_INV the first record will have matching source and target, for INV_dup it will be the second (indexing accordingly)
            self.assertTrue(
                records[int(sv_type == 'dup_INV')]['target_s'] == records[int(sv_type == 'dup_INV')]['source_s'])
            self.assertTrue(
                records[int(sv_type == 'dup_INV')]['target_e'] == records[int(sv_type == 'dup_INV')]['source_e'])
            if sv_type == 'dup_INV':
                self.assertTrue(records[0]['target_s'] == records[0]['target_e'] == '8')
                self.assertTrue(records[0]['ev_type'] == 'INVDUP')
                self.assertTrue(records[1]['ev_type'] == 'INV')
            else:
                self.assertTrue(records[1]['target_s'] == records[1]['target_e'] == '0')
                self.assertTrue(records[0]['ev_type'] == 'INV')
                self.assertTrue(records[1]['ev_type'] == 'INVDUP')

    def test_export_bedpe_INVdup(self):
        record = self.initialize_test(self.test_objects_INVdup, 'INVdup')[0]
        self.singleton_event_bed_tests([record], 'INVdup', 'chr19', '4')
        self.assertTrue((record['source_s'], record['source_e']) == (record['target_s'], record['target_e']) == ('0', '4'))
        self.assertTrue(record['ev_type'] == 'INVDUP')

if __name__ == "__main__":
    unittest.main()
