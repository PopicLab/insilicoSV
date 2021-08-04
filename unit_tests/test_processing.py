import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from processing import FormatterIO

class TestProcessing(unittest.TestCase):

    def setUp(self):
        # runs before every test
        self.formatter = FormatterIO("unit_tests/inputs/test.fna")

    def test_find_lcs(self):
        # if lowercase letters present in target, then take those letters in lcs
        self.assertEqual(self.formatter.find_lcs("ABC","cbC", input_is_list=False), list("bC"))

        # _ should always be included in lcs
        self.assertEqual(self.formatter.find_lcs(["A","_1","B"],["B","_1","A"]), ["_1"])
        self.assertEqual(self.formatter.find_lcs(["A","_1","B","C","D","E","F","A","_2"],["A","B","C","D","E","F","_1","A","_2"]), ["A","_1","A","_2"])
        with self.assertRaises(Exception):
            self.formatter.find_lcs("A_BC","ABF_A_", input_is_list=False)

        # insertions
        self.assertEqual(self.formatter.find_lcs(["A","_1","B"],["B","F","_1","A","C","D"]), ["_1"])
        
        # miscellaneous
        self.assertEqual(self.formatter.find_lcs([""],[""]), [])
    
    def test_export_bedpe(self):
        # Important: represents ground truth file
        pass

        

if __name__ == "__main__":
    unittest.main()