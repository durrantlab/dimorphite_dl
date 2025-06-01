import pytest

from dimorphite_dl import MoleculeRecord


@pytest.mark.parametrize(
    ("input_smiles", "exp_azides", "exp_neutral", "exp_canonical"),
    [
        ("C#CCO", "C#CCO", "C#CCO", "C#CCO"),
        ("Brc1cc[nH+]cc1", "Brc1cc[nH+]cc1", "Brc1ccncc1", "Brc1ccncc1"),
        ("C-N=[N+]=[N@H]", "C-N=[N+]=[N@H]", "CN=[N+]=N", "CN=[N+]=N"),
        ("O=P(O)(O)OCCCC", "O=P(O)(O)OCCCC", "CCCCOP(=O)(O)O", "CCCCOP(=O)(O)O"),
    ],
)
def test_molecule_preparation_steps(
    input_smiles, exp_azides, exp_neutral, exp_canonical
):
    mol = MoleculeRecord(input_smiles)
    mol.process_azides()
    assert mol.smiles == exp_azides, (
        f"after process_azides: got {mol.smiles!r}, expected {exp_azides!r}"
    )
    mol.neutralize()
    assert mol.smiles == exp_neutral, (
        f"after neutralize: got {mol.smiles!r}, expected {exp_neutral!r}"
    )
    mol.make_canonical()
    assert mol.smiles == exp_canonical, (
        f"after make_canonical: got {mol.smiles!r}, expected {exp_canonical!r}"
    )
