import pytest
import math


from insilicosv.variant_set import FromGrammarVariantSet 

# Dummy config to format error messages
DUMMY_VSET_CONFIG = {"config_descr": "test_sv_config"}

class TestResolverOrder:
    def test_independent(self):
        ranges = [[1, 2], [3, 4], [5, 6]]
        order = FromGrammarVariantSet.resolver_order(ranges, {"A": 0, "B": 1, "C": 2}, DUMMY_VSET_CONFIG)
        assert order == [0, 1, 2]
 
    def test_single_dependency(self):
        ranges = [["B", 5], [5, 5]]
        li = {"A": 0, "B": 1}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order == [1, 0]
 
    def test_transitive(self):
        ranges = [["2B", "2B"], ["C+1", "C+1"], [5, 5]]
        li = {"A": 0, "B": 1, "C": 2}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order == [2, 1, 0]
    
    def test_two_chains(self):
        ranges = [["2B", "2B"], ["C+1", "C+1"], [5, 5], ["C", 5], [4, "D"]]
        li = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order == [2, 1, 0, 3, 4]
 
    def test_double_dependency(self):
        ranges = [["B+C", 3], [2, 2], [3, 3]]
        li = {"A": 0, "B": 1, "C": 2}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order[-1] == 0

    def test_implicit_multiplication(self):
        ranges = [[3, 3], [4, 4], ["2AB", 5]]
        li = {"A": 0, "B": 1, "C": 2}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order[-1] == 2
 
    def test_min_and_max_different_letters(self):
        ranges = [[3, 3], ["A", "C"], [1, 1]]
        li = {"A": 0, "B": 1, "C": 2}
        order = FromGrammarVariantSet.resolver_order(ranges, li, DUMMY_VSET_CONFIG)
        assert order[-1] == 1
 
    def test_missing_symbol(self):
        with pytest.raises(ValueError):
            FromGrammarVariantSet.resolver_order(
                [[1, 1], ["Z", "Z"]], {"A": 0, "B": 1}, DUMMY_VSET_CONFIG
            )
 
    def test_cycle(self):
        with pytest.raises(SyntaxError):
            FromGrammarVariantSet.resolver_order(
                [["B", "B"], ["A", "A"]], {"A": 0, "B": 1}, DUMMY_VSET_CONFIG
            )
 

class TestEvalFormula:
    def test_non_string(self):
        assert FromGrammarVariantSet._eval_formula(7, {}, True, DUMMY_VSET_CONFIG) == 7
        assert FromGrammarVariantSet._eval_formula(None, {}, False, DUMMY_VSET_CONFIG) is None
 
    def test_letter(self):
        assert FromGrammarVariantSet._eval_formula("A", {"A": 6}, True, DUMMY_VSET_CONFIG) == 6
 
    def test_multiplication(self):
        assert FromGrammarVariantSet._eval_formula("2A", {"A": 6}, False, DUMMY_VSET_CONFIG) == 12
 
    def test_min_bound_division(self):
        assert FromGrammarVariantSet._eval_formula("A/3", {"A": 10}, True, DUMMY_VSET_CONFIG) == 4
 
    def test_max_bound_division(self):
        assert FromGrammarVariantSet._eval_formula("A/3", {"A": 10}, False, DUMMY_VSET_CONFIG) == 3
 
    def test_dependency_none(self):
        with pytest.raises(SyntaxError):
            FromGrammarVariantSet._eval_formula("A", {"A": None}, True, DUMMY_VSET_CONFIG)
 
    def test_letters_multiplication(self):
        assert FromGrammarVariantSet._eval_formula("2AB", {"A": 4, "B": 3}, False, DUMMY_VSET_CONFIG) == 24

    def test_letters_addition(self):
        assert FromGrammarVariantSet._eval_formula("A+ B", {"A": 4, "B": 3}, False, DUMMY_VSET_CONFIG) == 7

    def test_letters_substraction(self):
        assert FromGrammarVariantSet._eval_formula("A -B", {"A": 4, "B": 3}, False, DUMMY_VSET_CONFIG) == 1


class TestPickSymbolLengths:
    def test_fixed_ranges(self):
        ranges = [[5, 5], [7, 7]]
        order = [0, 1]
        num_letters = 2
        indexes_letter = {0: "A", 1: "B"}
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(ranges, order, num_letters, indexes_letter, DUMMY_VSET_CONFIG)
        assert lengths == [5, 7]
        assert min_lengths == [None, None]
 
    def test_transitive(self):
        ranges = [["2B", "2B"], ["C+1", "C+1"], [5, 5]]
        order = [2, 1, 0]
        num_letters = 3
        indexes_letter = {0: "A", 1: "B", 2: "C"}
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(ranges, order, num_letters, indexes_letter, DUMMY_VSET_CONFIG)

        assert lengths == [12, 6, 5]
        assert min_lengths == [None, None, None]
 
    def test_fraction(self):
        ranges = [[10, 10], ["A/3", 4]]
        order = [0, 1]
        num_letters = 2
        indexes_letter = {0: "A", 1: "B"}
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(ranges, order, num_letters, indexes_letter, DUMMY_VSET_CONFIG)

        assert lengths == [10, 4]
        assert min_lengths == [None, None]
 
    def test_dispersion(self):
        ranges = [[10, 10], [100, None]]
        order = [0, 1]
        num_letters = 1
        indexes_letter = {0: "A"}
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(ranges, order, num_letters, indexes_letter, DUMMY_VSET_CONFIG)

        assert lengths == [10, None]
        assert min_lengths == [None, 100]       
 
    def test_randomness(self):
        ranges = [[10, 20], ["A + 5", "2A"], ["B/2", "B"]]
        order = [0, 1, 2]
        num_letters = 3
        indexes_letter = {0: "A", 1: "B", 2: "C"}

        seen_lengths = set()
        for _ in range(100):
            lengths, _ = FromGrammarVariantSet.pick_symbol_lengths(ranges, order, num_letters, indexes_letter, DUMMY_VSET_CONFIG)
            seen_lengths.add(lengths[0])
            assert 10 <= lengths[0] <= 20
            assert lengths[0] + 5 <= lengths[1] <= 2 * lengths[0]
            assert math.ceil(lengths[1]/2) <= lengths[2] <= lengths[1]
        assert len(seen_lengths) > 1
        