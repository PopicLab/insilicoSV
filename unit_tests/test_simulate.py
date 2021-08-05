import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from simulate import SV_Simulator
from pysam import FastaFile
import yaml

class RandomSim():
    def __init__(self, increase, mode = None):
        # mode: None, "min," "mid," or "max" - specifies which number to pick in randint method
        # mode = None means to rely on self.curr_pos
        # replaces the random module in the simulator to make the results more predictable
        assert mode in [None, "min", "mid", "max"]

        self.curr_pos = 0
        self.increase = increase
        self.mode = mode
    def randint(self, min_int, max_int):
        assert (isinstance(min_int, int) and isinstance(max_int, int))
        if self.mode == None:
            if self.curr_pos + self.increase > max_int:
                return min_int
            tmp = self.curr_pos
            self.curr_pos += self.increase
            print("randint({}, {}) -> {}".format(min_int, max_int, tmp))
            return tmp
        elif self.mode == "min":
            return min_int
        elif self.mode == "max":
            return max_int
        elif self.mode == "mid":
            return int((min_int + max_int)/2)
        else:
            raise Exception ("Mode {} not valid",format(self.mode))

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
        def rm(file):
            if os.path.exists(file):
                os.remove(file)
            else:
                "File {} not detected".format(file)
        
        files = [self.ref, self.ref + ".fai", self.par, self.hap1, self.hap1 + ".fai", self.hap2, self.hap2 + ".fai", self.bed]  # remove reference's index (.fai) file because new one should be generated
        for file in files:
            rm(file)
    
    def get_actual_frag(self):
        fasta_out = FastaFile(self.hap1)  # also generates a .fai file
        return fasta_out.fetch(fasta_out.references[0], 0, fasta_out.get_reference_length(fasta_out.references[0]))

class TestSVSimulator(unittest.TestCase):

    def setUp(self):

        # runs before every test
        ref_file = "unit_tests/inputs/test.fna"
        par = "unit_tests/inputs/par.yaml"

        hap1 = "unit_tests/inputs/test1.fna"
        hap2 = "unit_tests/inputs/test2.fna"
        bed = "unit_tests/inputs/out.bed"

        self.test_objects_no_dis = [TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                        [par, {"SVs": [{"type": "delINVdup", "number": 1, "max_length": 5, "min_length": 5}]}],
                                        hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                                [par, {"SVs": [{"type": "delINVdup", "number": 1, "max_length": 5, "min_length": 5},
                                                                            {"type": "delINVdel", "number": 1, "min_length": 5, "max_length": 5},
                                                                            {"type": "dupINVdup", "number": 1, "min_length": 5, "max_length": 5}
                                                                            ]}],
                                                hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGAGTCAGGGAGCAAAAAAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                                [par, {"SVs": [{"type": "dupINVdel", "number": 1, "max_length": 5, "min_length": 5},
                                                                {"type": "delINV", "number": 1, "min_length": 5, "max_length": 5},
                                                                {"type": "INVdel", "number": 1, "min_length": 5, "max_length": 5}
                                                                ]}],
                                                hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                                [par, {"SVs": [{"type": "dupINV", "number": 1, "max_length": 5, "min_length": 5},
                                                                {"type": "INVdup", "number": 1, "min_length": 5, "max_length": 5}
                                                                ]}],
                                                hap1, hap2, bed)]
        self.test_objects_with_dis = [TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                        [par, {"SVs": [{"type": "TRA", "number": 2, "min_length": 5, "max_length": 5}]}],
                                        hap1, hap2, bed),
                                    TestObject([ref_file, {"Chromosome19": "CTCCGTCGTACTAGACAGGGTATATGTCTGTGTCTCAGTGAGACACTTAGCATGCAACTCAGTCTGTACTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT"}], 
                                        [par, {"SVs": [{"type": "dDUP-iDEL", "number": 1, "min_length": 5, "max_length": 5},
                                                        {"type": "INS-iDEL", "number": 1, "min_length": 5, "max_length": 5},
                                                        {"type": "dDUP", "number": 1, "min_length": 5, "max_length": 5}]}],
                                        hap1, hap2, bed)]



    def test_produce_variant_genome(self):

        # ====== Test SVs without dispersions ======
        # also tests whether overlapping detected

        # dupINVdup, delINVdel, delINVdup
        config_no_dis = self.test_objects_no_dis[1]
        config_no_dis.initialize_files()    # generate reference and parameter files
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par, random_gen = RandomSim(0, "max"))   # random_gen ensures that SVs are all homozygous
        self.assertEqual(curr_sim.produce_variant_genome(config_no_dis.hap1, config_no_dis.hap2, config_no_dis.bed, random_gen = RandomSim(10)), True)
        changed_frag = config_no_dis.get_actual_frag()   # fetch fragment produced by simulator to compare with answer
        config_no_dis.remove_test_files()   # remove any output files and .fai files
        # Work
        # CTCCG TCGTA CTAGA CAGCT CCCGA CAGAG CACTG GTGTC TTGTT TCTTT AAACA CCAGT ATTTA GATGCACTATCTCTCCGT
        # TCTAG TACGA CTAGA CAGCT _____ CTCTG _____ GTGTC TTGTT TGTTTAAAGAAACAA AAACA CCAGT ATTTA GATGCACTATCTCTCCGT
        self.assertEqual(changed_frag, "TCTAGTACGACTAGACAGCTCTCTGGTGTCTTGTTTGTTTAAAGAAACAAAAACACCAGTATTTAGATGCACTATCTCTCCGT")
        
        # dupINVdel, delINV, INVdel
        config_no_dis = self.test_objects_no_dis[2]
        config_no_dis.initialize_files()
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par, random_gen = RandomSim(0, "max"))
        self.assertEqual(curr_sim.produce_variant_genome(config_no_dis.hap1, config_no_dis.hap2, config_no_dis.bed, random_gen = RandomSim(10)), True)
        changed_frag = config_no_dis.get_actual_frag()
        config_no_dis.remove_test_files()
        # CTCCG TCGTA CTAGA CAGCT CCCGA GTCAG GGAGC AAAAA AGTGT GACAC TAGTC CACAG GTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        # CTCCG TACGA CGGAG CAGCT _____ CTGAC GCTCC _____ AGTGT GACAC TAGTC CACAG GTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        self.assertEqual(changed_frag, "CTCCGTACGACGGAGCAGCTCTGACGCTCCAGTGTGACACTAGTCCACAGGTGAGAAACACAAATATTCAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT")
        #ACACTAGTCCACAGGTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        
        # dupINV, INVdup
        config_no_dis = self.test_objects_no_dis[3]
        config_no_dis.initialize_files()
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par, random_gen = RandomSim(0, "max"))
        self.assertEqual(curr_sim.produce_variant_genome(config_no_dis.hap1, config_no_dis.hap2, config_no_dis.bed, random_gen = RandomSim(10)), True)
        changed_frag = config_no_dis.get_actual_frag()
        config_no_dis.remove_test_files()
        # ACACT AGTCC ACAGG TGAGA ATCTT GTTTC TTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        # ACACT GGACTAGTGT TCTCACCTGT TGAGA ATCTT GTTTC TTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        self.assertEqual(changed_frag, "ACACTGGACTAGTGTTCTCACCTGTTGAGAATCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT")

        # ====== Test SVs with dispersions ======
        # TRA
        config_with_dis = self.test_objects_with_dis[0]
        config_with_dis.initialize_files()
        curr_sim = SV_Simulator(config_with_dis.ref, config_with_dis.par, random_gen = RandomSim(10, "max"))
        self.assertEqual(curr_sim.produce_variant_genome(config_with_dis.hap1, config_with_dis.hap2, config_with_dis.bed, random_gen = RandomSim(10)), True)
        changed_frag = config_with_dis.get_actual_frag()
        config_with_dis.remove_test_files()
        # CTCCG TCGTA CTAGA CAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        # TCTAG TACGA CTAGA CAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        self.assertEqual(changed_frag, "CTAGATCGTACTCCGCAGCTCACTGCAGAGCCCGAGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT")

        # dDUP-iDEL, INS-IDEL, dDUP
        config_with_dis = self.test_objects_with_dis[1]
        config_with_dis.initialize_files()
        curr_sim = SV_Simulator(config_with_dis.ref, config_with_dis.par, random_gen = RandomSim(10, "max"))
        self.assertEqual(curr_sim.produce_variant_genome(config_with_dis.hap1, config_with_dis.hap2, config_with_dis.bed, random_gen = RandomSim(5)), True)
        changed_frag = config_with_dis.get_actual_frag()
        config_with_dis.remove_test_files()
        # CTCCG TCGTA CTAGA CAGGG TATAT GTCTG TGTCT CAGTG AGACA CTTAGCATGCAACTCAGTCTGTACTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        # CTCCG _____ CTCCG TCGTA TATAT GTCTG TATATTGTCT CAGTG AGACA CTTAGCATGCAACTCAGTCTGTACTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        self.assertEqual(changed_frag, "CTCCGCTCCGTCGTATATATGTCTGTATATTGTCTCAGTGAGACACTTAGCATGCAACTCAGTCTGTACTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT")
    
    
    def test_choose_rand_pos(self):
        # test SVs without dispersions
        config_no_dis = self.test_objects_no_dis[0]
        config_no_dis.initialize_files()
        curr_sim = SV_Simulator(config_no_dis.ref, config_no_dis.par, random_gen = RandomSim(10, "max"))
        random_gen = RandomSim(10)
        self.assertEqual(curr_sim.choose_rand_pos(curr_sim.svs, curr_sim.ref_fasta, random_gen), {"Chromosome19": [(0,5), (5,10), (10,15)]})

if __name__ == "__main__":
    unittest.main()