import pytest

from insilicosv.variant_set_makers import Symbol

def test_symbol():
    assert Symbol('A') == Symbol('A')
    with pytest.raises(Exception):
        Symbol('a')
    with pytest.raises(Exception):
        Symbol("A'")
    with pytest.raises(Exception):
        Symbol("A*")
    with pytest.raises(Exception):
        Symbol("_")

    assert (sorted([Symbol('B'), Symbol('_1'), Symbol('A')]) ==
            [Symbol('A'), Symbol('B'), Symbol('_1')])
# end: def test_symbol()

