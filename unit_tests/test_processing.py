import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from processing import FormatterIO

class TestProcessing(unittest.TestCase):

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
        self.formatter = FormatterIO("unit_tests/inputs/test.fna")

    def tearDown(self):
        # runs after every test
        pass

    def test1_yaml(self):
        # primarily checks error detection system
        pass

    def test_find_lcs(self):
        # if lowercase letters present in target, then take those letters in lcs
        self.assertEqual(self.formatter.find_lcs("ABC","cbC"), "bC")

        # _ should always be included in lcs
        self.assertEqual(self.formatter.find_lcs("A_B","B_A"), "_")
        self.assertEqual(self.formatter.find_lcs("A_BCDEFA_","ABCDEF_A_"), "A_A_")
        with self.assertRaises(Exception):
            self.formatter.find_lcs("A_BC","ABF_A_")

        # insertions
        self.assertEqual(self.formatter.find_lcs("A_B","BF_ACD"), "_")
        
        # miscellaneous
        self.assertEqual(self.formatter.find_lcs("",""), "")

        

if __name__ == "__main__":
    unittest.main()