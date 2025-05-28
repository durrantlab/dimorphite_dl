import pytest
from rdkit import Chem

from dimorphite_dl import MoleculeRecord
from dimorphite_dl.protonate.detect import ProtonationSiteDetector


def compare_smiles(smiles1, smiles2):
    detected_can = Chem.MolToSmiles(Chem.MolFromSmiles(smiles1))
    expected_can = Chem.MolToSmiles(Chem.MolFromSmiles(smiles2))
    assert detected_can == expected_can, f"got {smiles1}, expected {smiles2}"


def compare_smarts(smarts1, smarts2):
    detected_can = Chem.MolToSmarts(Chem.MolFromSmarts(smarts1))
    expected_can = Chem.MolToSmarts(Chem.MolFromSmarts(smarts2))
    assert detected_can == expected_can, f"got {smarts1}, expected {smarts2}"


@pytest.mark.parametrize(
    ("smiles", "smiles_prepped_correct", "expected_smarts", "expected_idx"),
    [
        ("C#CCO", "[H]C#CC([H])([H])O[H]", "[C:1]-[O:2]-[#1]", (2, 3, 7)),
        ("Brc1cc[nH+]cc1", "[H]c1nc([H])c([H])c(Br)c1[H]", "[n&+0&H0:1]", (4,)),
        (
            "C-N=[N+]=[N@H]",
            "[H]N=[N+]=NC([H])([H])[H]",
            "[N&+0:1]=[N&+:2]=[N&+0:3]-[#1]",
            (1, 2, 3, 7),
        ),
        (
            "O=P(O)(O)OCCCC",
            "[H]OP(=O)(O[H])OC([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]",
            "[P&X4:1](=[O:2])(-[O&X2:3]-[#1])(-[O&+0:4])-[O&X2:5]-[#1]",
            (5, 6, 7, 18, 4, 8, 19),
        ),
    ],
)
def test_substructure_detect(
    smiles, smiles_prepped_correct, expected_smarts, expected_idx
):
    mol_record = MoleculeRecord(smiles)

    detector = ProtonationSiteDetector()

    # prepare molecule
    mol = mol_record.prepare_for_protonation()
    smiles_prepped = Chem.MolToSmiles(mol)
    compare_smiles(smiles_prepped, smiles_prepped_correct)

    # detect substructures
    substructures = detector._detect_all_sites_in_molecule(mol)
    sub_match = substructures[0]

    # instead of raw string equality, canonicalize both SMARTS and compare
    compare_smarts(sub_match.substructure.smarts, expected_smarts)
    # atom indices should still be the same
    assert sub_match.idx_atom == expected_idx
