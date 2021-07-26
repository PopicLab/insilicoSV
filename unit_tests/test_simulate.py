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
        self.sv = Structural_Variant(8, (10,10))

    def tearDown(self):
        # runs after every test
        pass

    def test_change_fragment(self):
        pass

        

if __name__ == "__main__":
    unittest.main()