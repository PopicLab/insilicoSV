import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from simulate import Structural_Variant

class TestStructuralVariant(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # runs before all tests
        pass

    @classmethod
    def tearDownClass(cls):
        # runs after all tests
        pass

    def setUp(self):
        # runs before every test
        self.sv_no_dispersion = Structural_Variant(8, [(5,5)])
        self.sv_with_dispersion = Structural_Variant(5, [(5,5)])

    def tearDown(self):
        # runs after every test
        pass

    def test_generate_blocks(self):
        self.assertEqual(self.sv_no_dispersion.generate_blocks(),[["c","b","C"]])

    def test_change_fragment(self):
        def initialize_event(event,start_pos, chr_id, frag):
            event.start = start_pos
            event.end = event.start + event.length
            event.source_chr = chr_id
            if not event.symbol.startswith("_"):
                event.source_frag = frag

        # non-dispersion events
        bases = ["A", "T", "G", "C"]
        base_idx = 0
        start_pos = 10
        for event in self.sv_no_dispersion.source_events:
            initialize_event(event, start_pos, "chr21", bases[base_idx] * event.length)
            start_pos += event.length
            base_idx += 1
        #ABC -> cbC
        #AAAAATTTTTGGGGG -> CCCCCAAAAAGGGGG
        self.assertEqual(self.sv_no_dispersion.change_fragment(), [["chr21", 10, 25, "CCCCCAAAAAGGGGG"]])

        # SVs with dispersion events
        base_idx = 0
        start_pos = 25
        for event in self.sv_with_dispersion.source_events:
            initialize_event(event, start_pos, "chr21", bases[base_idx] * event.length)
            start_pos += event.length
            base_idx += 1
        #A_B -> B_A
        #AAAAA_GGGGG -> GGGGG_AAAAA
        self.assertEqual(self.sv_with_dispersion.change_fragment(), [["chr21", 25, 30, "GGGGG"], ["chr21", 35, 40, "AAAAA"]])


class TestSVSimulator(unittest.TestCase):
    def test_produce_variant_genome(self):
        # test insertions

        # test SVs without dispersions

        # test SVs with dispersions
        pass

        

if __name__ == "__main__":
    unittest.main()