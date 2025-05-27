"""
Tests for loading protonation site data
"""

from dimorphite_dl.protonate import SubstructureData


def test_data_init():
    pka_data = SubstructureData()
    pka_data2 = SubstructureData()
    assert pka_data == pka_data2

    n_substructures = len(pka_data._data)
    assert n_substructures == 41
