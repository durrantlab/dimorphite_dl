from dimorphite_dl.protonate.data import PKaData


def test_data_init():
    pka_data = PKaData()
    pka_data2 = PKaData()
    assert pka_data == pka_data2

    n_substructures = len(pka_data._data)
    assert n_substructures == 41
