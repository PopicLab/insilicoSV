import pytest

from insilicosv.variant_set import Symbol, FromGrammarVariantSet, Operation, TransformType

def get_operations(operations):
    return set((ope.transform.transform_type, ope.transform.is_in_place, (ope.source_breakend_region.start_breakend if ope.source_breakend_region is not None else None,
                                                                          ope.source_breakend_region.end_breakend if ope.source_breakend_region is not None else None,
                                                                          ope.target_insertion_breakend,
                                                                          ope.target_insertion_order[1] if ope.target_insertion_order is not None else None))
                  for ope in operations)

def test_pick_symbol_length():
    letter_indexes = {'A': 0, 'B': 1, 'C': 2}
    length_ranges = [[5, 10], ['C*0.3', 'A-2'], ['A', 'A']]
    dispersion_ranges = [['B/10', 'B']]
    lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                     dispersion_ranges=dispersion_ranges,
                                                                     letter_indexes=letter_indexes)
    assert 5 <= lengths[0] <= 10
    assert lengths[0] - 2 >= lengths[1]
    assert lengths[2] * 0.3 <= lengths[1]
    assert lengths[2] == lengths[0]
    assert lengths[1]/10 <= lengths[3] <= lengths[1]

    length_ranges = [[50, 100], ['C*(0.3*A + 30)', '(A-2)*1000'], ['A', 'A']]
    dispersion_ranges = [['B/10', 'B']]
    lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                     dispersion_ranges=dispersion_ranges,
                                                                     letter_indexes=letter_indexes)
    assert 50 <= lengths[0] <= 100
    assert (lengths[0] - 2)*1000 >= lengths[1]
    assert lengths[2] * (0.3*lengths[0]+30) <= lengths[1]
    assert lengths[2] == lengths[0]
    assert lengths[1] / 10 <= lengths[3] <= lengths[1]

    valid = False
    try:
        length_ranges = [['B', 10], ['C*0.3', 'A-2'], ['A', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid

    valid = False
    try:
        length_ranges = [[11, 10], ['C*0.3', 'A-2'], ['A', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid

    length_ranges = [[5, 10], ['C*0.03A', 'A-2'], ['A', 'A']]
    dispersion_ranges = [['B/10', 'B']]
    lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                     dispersion_ranges=dispersion_ranges,
                                                                     letter_indexes=letter_indexes)
    assert 5 <= lengths[0] <= 10
    assert lengths[0] - 2 >= lengths[1]
    assert lengths[2] * 0.03 * lengths[0] <= lengths[1]
    assert lengths[2] == lengths[0]
    assert lengths[1] / 10 <= lengths[3] <= lengths[1]

    valid = False
    try:
        length_ranges = [[5, 10], ['(C*0.3', 'A-2'], ['A', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid

    valid = False
    try:
        length_ranges = [[5, 10], ['C(*0.3)', 'A-2'], ['A', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
        print('lengths', lengths)
    except:
        valid = True
    assert valid
    valid = False
    try:
        length_ranges = [[5, 10], ['C*0.3', 'A-2'], ['G', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid

    valid = False
    try:
        length_ranges = [[5, 10], ['C*0.3', 'A-2'], ['_', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid

    valid = False
    try:
        length_ranges = [[5, 10], ['A+50', 'A-2'], ['_', 'A']]
        dispersion_ranges = [['B/10', 'B']]
        lengths, min_lengths = FromGrammarVariantSet.pick_symbol_lengths(length_ranges=length_ranges,
                                                                         dispersion_ranges=dispersion_ranges,
                                                                         letter_indexes=letter_indexes)
    except:
        valid = True
    assert valid


def test_FromGrammarVariantSetMaker():
    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source':'ABC_', 'target': 'A_BC', 'length_ranges':[[1,1],
                                           [2,2], [3,3]], 'dispersion_ranges':[[1,1]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 2, 3, 1]
    operations = get_operations(sv.operations)
    print('operations', operations)
    assert set([(TransformType.IDENTITY, False, (1, 2, 4, 1)), (TransformType.IDENTITY, False, (2, 3, 4, 2)),
                (TransformType.DEL, True, (1, 2, None, None)), (TransformType.DEL, True, (2, 3, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'AB_C', 'target': 'A_BC', 'length_ranges': [[1, 1],
                                                                                                              [2, 2],
                                                                                                              [3, 3]],
                                                        'dispersion_ranges': [[1, 1]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 2, 1, 3]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (1, 2, 3, 1)), (TransformType.DEL, True, (1, 2, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'A_B_C', 'target': 'c_A_B', 'length_ranges': [[1, 1],
                                                                                                              [2, 2],
                                                                                                              [3, 3]],
                                                        'dispersion_ranges': [[1, 1], [1, 1]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 1, 2, 1, 3]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (0, 1, 2, 2)), (TransformType.IDENTITY, False, (2, 3, 4, 3)),
                (TransformType.INV, False, (4, 5, 0, 1)),
                (TransformType.DEL, True, (0, 1, None, None)), (TransformType.DEL, True, (2, 3, None, None)),
                (TransformType.DEL, True, (4, 5, None, None))
                ]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'nrTRA', 'number': 1,
                                                        'length_ranges': [[1, 1]],
                                                        'dispersion_ranges': [[2, 2]]})
    sv = set_maker.simulate_sv()
    assert (sv.breakend_interval_lengths == [1, 2]) or (sv.breakend_interval_lengths == [2, 1])
    operations = get_operations(sv.operations)
    assert ((set([(TransformType.IDENTITY, False, (0, 1, 2, 1)), (TransformType.DEL, True, (0, 1, None, None))]) == operations) or
            (set([(TransformType.IDENTITY, False, (1, 2, 0, 1)), (TransformType.DEL, True, (1, 2, None, None))]) == operations))

    set_maker = FromGrammarVariantSet(vset_config={'type': 'rTRA', 'number': 1,
                                                        'length_ranges': [[1, 1], [3,3]],
                                                        'dispersion_ranges': [[2, 2]]})
    sv = set_maker.simulate_sv()
    assert (sv.breakend_interval_lengths == [1, 2, 3]) or (sv.breakend_interval_lengths == [3, 2, 1])
    operations = get_operations(sv.operations)
    assert (set([(TransformType.IDENTITY, False, (0, 1, 2, 2)), (TransformType.DEL, True, (0, 1, None, None)),
                 (TransformType.IDENTITY, False, (2, 3, 0, 1)), (TransformType.DEL, True, (2, 3, None, None))]) == operations)

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'A___B', 'target': '_A_B_', 'length_ranges': [[1, 1],
                                                                                                              [3, 3]],
                                                        'dispersion_ranges': [[1, 1], [5, 5], [3, 3]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 1, 5, 3, 3]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (0, 1, 2, 1)), (TransformType.IDENTITY, False, (4, 5, 3, 2)),
                (TransformType.DEL, True, (0, 1, None, None)), (TransformType.DEL, True, (4, 5, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'AB', 'target': 'ACB', 'length_ranges': [[1, 1], [2, 2],
                                                                                                                [3, 3]],
                                                        'dispersion_ranges': []})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 2]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (None, None, 1, 1))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'ABC', 'target': 'BCA', 'length_ranges': [[1, 1],
                                                                                                              [2, 2],
                                                                                                              [3, 3]],
                                                        'dispersion_ranges': []})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 2, 3]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (0, 1, 3, 1)), (TransformType.DEL, True, (0, 1, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'AB_', 'target': '_bA', 'length_ranges': [[1, 1],
                                                                                                            [3, 3]],
                                                        'dispersion_ranges': [[4, 4]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 3, 4]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, False, (0, 1, 3, 2)), (TransformType.INV, False, (1, 2, 3, 1)),
                (TransformType.DEL, True, (0, 1, None, None)), (TransformType.DEL, True, (1, 2, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'AB', 'target': 'A',
                                                        'length_ranges': [[1, 1], [2, 2]],
                                                        'dispersion_ranges': []})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 2]
    operations = get_operations(sv.operations)
    assert set([(TransformType.DEL, True, (1, 2, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'Custom', 'number': 1,
                                                        'source': 'A___B', 'target': '__B_C', 'length_ranges': [[1, 1],
                                                                                                                [3, 3],
                                                                                                                [4, 4]],
                                                        'dispersion_ranges': [[1, 1], [5, 5], [3, 3]]})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1, 1, 5, 3, 3]
    operations = get_operations(sv.operations)
    assert set(
        [(TransformType.DEL, True, (0, 1, None, None)), (TransformType.IDENTITY, False, (4, 5, 3, 1)),
         (TransformType.DEL, True, (4, 5, None, None)),
         (TransformType.IDENTITY, False, (None, None, 4, 2))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'DEL', 'length_ranges': [[3, 3]], 'number': 1})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [3]
    operations = get_operations(sv.operations)
    assert set(
        [(TransformType.DEL, True, (0, 1, None, None))]) == operations

    set_maker = FromGrammarVariantSet(vset_config={'type': 'SNP', 'number': 1})
    sv = set_maker.simulate_sv()
    assert sv.breakend_interval_lengths == [1]
    operations = get_operations(sv.operations)
    assert set([(TransformType.IDENTITY, True, (0, 1, None, None))]) == operations

if __name__ == '__main__':
    test_pick_symbol_length()