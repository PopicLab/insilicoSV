import copy
import math
import pytest


from insilicosv.utils import pick_symbol_lengths 

# Dummy config to format error messages
DUMMY_VSET_CONFIG = {"config_descr": "test_sv_config"}

class TestPickSymbolLengths:

    def test_complex_dependencies_variance_and_isolation(self):
        """
        Tests that when simulating multiple SVs:
        - The generated lengths are correct.
        - Randomness works.
        - The bounds are not modified from one iteration to the next.
        """
        original_ranges = [[10, 20], ["A + 5", "2A"], ["B/2", "B"]]
        dispersion_ranges = []
        letter_indexes = {"A": 0, "B": 1, "C": 2}
        
        input_ranges = copy.deepcopy(original_ranges)

        results_A = set()
        results_B = set()

        for _ in range(100):
            lengths, min_lengths = pick_symbol_lengths(
                input_ranges, dispersion_ranges, letter_indexes, DUMMY_VSET_CONFIG
            )
            
            assert input_ranges == original_ranges, "Function mutated the original length_ranges!"

            A, B, C = lengths
            
            assert 10 <= A <= 20
            assert A + 5 <= B <= 2 * A
            assert math.ceil(B / 2) <= C <= math.floor(B)

            results_A.add(A)
            results_B.add(B)

        assert len(results_A) > 1, "Randomness failed: Symbol A was identical across 100 iterations."
        assert len(results_B) > 1, "Randomness failed: Symbol B was identical across 100 iterations."

    def test_transitive_and_implicit_math(self):
        """
        Tests that the deque processes the order of the know variable doesn't matter.
        Test implicit math notations
        """
        length_ranges = [["2B", "2B"], ["C+1", "C+1"], [5, 5]]
        letter_indexes = {"A": 0, "B": 1, "C": 2}

        lengths, _ = pick_symbol_lengths(
            length_ranges, [], letter_indexes, DUMMY_VSET_CONFIG
        )

        assert lengths == [12, 6, 5]

    def test_fractional_bounds_rounding(self):
        """
        Check the rounding of math expressions
        """
        length_ranges = [[10, 10], ["A/3", "A/2"]] 
        letter_indexes = {"A": 0, "B": 1}

        lengths, _ = pick_symbol_lengths(
            length_ranges, [], letter_indexes, DUMMY_VSET_CONFIG
        )

        assert lengths[0] == 10
        assert lengths[1] in [4, 5]

    @pytest.mark.parametrize("length_ranges, letter_indexes, expected_error, match_text", [
        #Check failure modes.

        # Missing bounds
        ([[10]], {"A": 0}, SyntaxError, r"must be a list of \[min, max\] pairs"),
        
        # Min > Max
        ([[20, 10]], {"A": 0}, SyntaxError, "max bound less than min bound"),
        
        # Min > Max from math expression
        ([[10, 10], ["A", "A/2"]], {"A": 0, "B": 1}, SyntaxError, "max bound less than min bound"),
        
        # Negative length
        ([[-5, 5]], {"A": 0}, ValueError, "min length cannot be negative"),
        
        # Negative length from math expression
        ([[10, 10], ["A-20", "A"]], {"A": 0, "B": 1}, SyntaxError, "yielded negative bound"),
        
        # Missing Dependency
        ([[10, 10], ["Z", "Z"]], {"A": 0, "B": 1}, ValueError, "one of which is not define"),
        
        # Cyclic Dependency
        ([["B", "B"], ["A", "A"]], {"A": 0, "B": 1}, SyntaxError, "cyclic dependency"),
        
        # Invalid Math
        ([[10, 10], ["A/0", "A"]], {"A": 0, "B": 1}, SyntaxError, "Invalid math operation"),
        
        # Dependency on an Unbounded Symbol
        ([[None, None], ["A", "A"]], {"A": 0, "B": 1}, SyntaxError, "cannot depend on an unbounded symbol"),

    ], ids=[
        "malformed_pair", "hardcoded_min_max", "evaluated_min_max", "hardcoded_negative",
        "evaluated_negative", "undefined_symbol", "cyclic_dependency", "math_error", "unbounded_dependency"
    ])
    
    def test_pick_symbol_lengths_failure_modes(self, length_ranges, letter_indexes, expected_error, match_text):
        """Tests that all invalid configurations safely crash."""
        
        with pytest.raises(expected_error, match=match_text):
            pick_symbol_lengths(
                length_ranges, 
                dispersion_ranges=[], 
                letter_indexes=letter_indexes, 
                vset_config=DUMMY_VSET_CONFIG
            )