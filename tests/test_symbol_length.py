import copy
import math
import pytest

# Adjust this import based on your exact project structure. 
# It targets the class where pick_symbol_lengths is defined.
from insilicosv.variant_set import FromGrammarVariantSet 

# A dummy config required by the function to format error messages
DUMMY_VSET_CONFIG = {"config_descr": "test_sv_config"}

class TestPickSymbolLengths:

    def test_complex_dependencies_variance_and_isolation(self):
        """
        Tests that when called multiple times (simulating multiple SVs):
        1. The generated lengths strictly adhere to the mathematical dependencies.
        2. The values VARY between iterations (randomness works).
        3. The original configuration object is NEVER mutated (state isolation).
        """
        # A: [10, 20]
        # B: [A + 5, 2 * A]
        # C: [B / 2, B]
        original_ranges = [[10, 20], ["A + 5", "2A"], ["B/2", "B"]]
        dispersion_ranges = []
        letter_indexes = {"A": 0, "B": 1, "C": 2}
        
        # Deepcopy to ensure we can verify the function doesn't mutate inputs
        input_ranges = copy.deepcopy(original_ranges)

        results_A = set()
        results_B = set()

        for _ in range(100): # Simulate generating 100 SVs
            lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(
                input_ranges, dispersion_ranges, letter_indexes, DUMMY_VSET_CONFIG
            )
            
            # 1. Assert original state was NOT mutated
            assert input_ranges == original_ranges, "Function mutated the original length_ranges!"

            # 2. Assert constraints are perfectly maintained
            A, B, C = lengths
            
            assert 10 <= A <= 20
            assert A + 5 <= B <= 2 * A
            
            # Fractional limits use ceil for min and floor for max
            assert math.ceil(B / 2) <= C <= math.floor(B)

            # Store for variance check
            results_A.add(A)
            results_B.add(B)

        # 3. Assert variance (Values shouldn't be identical across 100 runs)
        assert len(results_A) > 1, "Randomness failed: Symbol A was identical across 100 iterations."
        assert len(results_B) > 1, "Randomness failed: Symbol B was identical across 100 iterations."

    def test_transitive_and_implicit_math(self):
        """
        Tests that the deque processes out-of-order transitive dependencies 
        (A depends on B, B depends on C, C is known) and parses implicit multiplication (2A = 2*A).
        """
        # A depends on 2B. B depends on C + 1. C is exactly 5.
        length_ranges = [["2B", "2B"], ["C+1", "C+1"], [5, 5]]
        letter_indexes = {"A": 0, "B": 1, "C": 2}

        lengths, _ = FromGrammarVariantSet.pick_symbol_lengths(
            length_ranges, [], letter_indexes, DUMMY_VSET_CONFIG
        )

        assert lengths == [12, 6, 5] # C=5, B=(5+1)=6, A=(2*6)=12

    def test_fractional_bounds_rounding(self):
        """
        Ensures that fractional math expressions round correctly:
        Min bounds should round UP (ceil), Max bounds should round DOWN (floor)
        to ensure the integer stays strictly inside the fraction constraint.
        """
        # 10/3 = 3.33 -> ceil -> 4
        # 10/2 = 5.0  -> floor -> 5
        length_ranges = [[10, 10], ["A/3", "A/2"]] 
        letter_indexes = {"A": 0, "B": 1}

        lengths, _ = FromGrammarVariantSet.pick_symbol_lengths(
            length_ranges, [], letter_indexes, DUMMY_VSET_CONFIG
        )

        assert lengths[0] == 10
        assert lengths[1] in [4, 5]


    # ==========================================
    # EXHAUSTIVE FAILURE MODE TESTING
    # ==========================================

    @pytest.mark.parametrize("length_ranges, letter_indexes, expected_error, match_text", [
        
        # 1. Malformed bounds format (List of 1 instead of 2)
        # CHANGED: Expected SyntaxError instead of ValueError, and used raw string r"..."
        ([[10]], {"A": 0}, SyntaxError, r"must be a list of \[min, max\] pairs"),
        
        # 2. Hardcoded Min > Max
        ([[20, 10]], {"A": 0}, SyntaxError, "max bound less than min bound"),
        
        # 3. Evaluated Min > Max (e.g., A=10 -> B in [10, 5])
        ([[10, 10], ["A", "A/2"]], {"A": 0, "B": 1}, SyntaxError, "max bound less than min bound"),
        
        # 4. Hardcoded Negative length
        ([[-5, 5]], {"A": 0}, ValueError, "min length cannot be negative"),
        
        # 5. Evaluated Negative Length (A=10, B = A-20 = -10)
        # CHANGED: Expected SyntaxError instead of ValueError because the code wraps it in a try/except
        ([[10, 10], ["A-20", "A"]], {"A": 0, "B": 1}, SyntaxError, "yielded negative bound"),
        
        # 6. Missing Dependency / Undefined Symbol ("Z" doesn't exist)
        ([[10, 10], ["Z", "Z"]], {"A": 0, "B": 1}, ValueError, "one of which is not define"),
        
        # 7. Cyclic Dependency (A depends on B, B depends on A)
        ([["B", "B"], ["A", "A"]], {"A": 0, "B": 1}, SyntaxError, "cyclic dependency"),
        
        # 8. Invalid Math execution (Division by Zero)
        ([[10, 10], ["A/0", "A"]], {"A": 0, "B": 1}, SyntaxError, "Invalid math operation"),
        
        # 9. Dependency on an Unbounded Symbol (A is [None, None])
        ([[None, None], ["A", "A"]], {"A": 0, "B": 1}, SyntaxError, "cannot depend on an unbounded symbol"),

    ], ids=[
        "malformed_pair", "hardcoded_min_max", "evaluated_min_max", "hardcoded_negative",
        "evaluated_negative", "undefined_symbol", "cyclic_dependency", "math_error", "unbounded_dependency"
    ])
    
    def test_pick_symbol_lengths_failure_modes(self, length_ranges, letter_indexes, expected_error, match_text):
        """Tests that all invalid configurations safely crash with the exact correct insilicoSV error."""
        
        with pytest.raises(expected_error, match=match_text):
            FromGrammarVariantSet.pick_symbol_lengths(
                length_ranges, 
                dispersion_ranges=[], 
                letter_indexes=letter_indexes, 
                vset_config=DUMMY_VSET_CONFIG
            )