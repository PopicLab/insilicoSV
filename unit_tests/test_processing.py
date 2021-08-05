import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from processing import FormatterIO

class TestProcessing(unittest.TestCase):

    def setUp(self):
        # runs before every test
        self.formatter = FormatterIO("unit_tests/inputs/par.yaml")

    def test_export_bedpe(self):
        # Important: represents ground truth file
        pass

        

if __name__ == "__main__":
    unittest.main()