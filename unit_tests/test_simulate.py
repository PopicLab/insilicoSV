import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from simulate import SV_Simulator
from pysam import FastaFile

class RandomSim():
    def __init__(self, increase, mode = None):
        # mode: None, "min," "mid," or "max" - specifies which number to pick in randint method
        # mode = None means to rely on self.curr_pos
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
            return tmp
        elif self.mode == "min":
            return min_int
        elif self.mode == "max":
            return max_int
        elif self.mode == "mid":
            return int((min_int + max_int)/2)
        else:
            raise Exception ("Mode {} not valid",format(self.mode))


class TestSVSimulator(unittest.TestCase):

    def setUp(self):
        # runs before every test
        self.ref_file1 = "unit_tests/inputs/test.fna"
        self.par1 = "unit_tests/inputs/par1.yaml"
        self.fna_1_1 = "unit_tests/inputs/test_1.fna"
        self.fna_2_1 = "unit_tests/inputs/test_2.fna"
        self.bed1 = "unit_tests/inputs/out.bed"
        self.sim1 = SV_Simulator(self.ref_file1, self.par1, random_gen = RandomSim(10))
    
    def tearDown(self):
        #runs after every test
        self.sim1.close()

    def test_produce_variant_genome(self):

        # test SVs without dispersions
        curr_sim = self.sim1
        random_gen = RandomSim(10)
        self.assertEqual(curr_sim.produce_variant_genome(self.fna_1_1, self.fna_2_1, self.bed1, random_gen = random_gen), True)
        fasta_out1 = FastaFile(self.fna_1_1)
        changed_frag = fasta_out1.fetch(fasta_out1.references[0], 0, fasta_out1.get_reference_length(fasta_out1.references[0]))
        # CTCCG TCGTA CTAGA CAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT
        # TCTAG TACGA CTAGA CAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT

        self.assertEqual(changed_frag, "TCTAGTACGACTAGACAGCTCCCGACAGAGCACTGGTGTCTTGTTTCTTTAAACACCAGTATTTAGATGCACTATCTCTCCGT")

        # test SVs with dispersions

        # test insertions
    
    def test_rand_select_svs(self):
        # test SVs without dispersions
        curr_sim = self.sim1
        random_gen = RandomSim(10)
        self.assertEqual(curr_sim.rand_select_svs(curr_sim.svs, curr_sim.ref_fasta, random_gen), [(0,5), (5,10), (10,15)])


        

if __name__ == "__main__":
    unittest.main()