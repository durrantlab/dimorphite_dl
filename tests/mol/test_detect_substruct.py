import pytest
from conftest import compare_smarts, compare_smiles  # type: ignore
from rdkit import Chem

from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.protonate.detect import ProtonationSiteDetector
from dimorphite_dl.protonate.site import ProtonationSite


@pytest.mark.parametrize(
    ("smiles", "smiles_prepped_correct", "expected_smarts", "expected_idxs_match"),
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
    smiles, smiles_prepped_correct, expected_smarts, expected_idxs_match
):
    mol_record = MoleculeRecord(smiles)

    detector = ProtonationSiteDetector()

    # prepare molecule
    mol = mol_record.prepare_for_protonation()
    smiles_prepped = Chem.MolToSmiles(mol)
    compare_smiles(smiles_prepped, smiles_prepped_correct)

    # detect substructures
    substructures = list(detector._detect_all_sites_in_molecule(mol))
    sub_match = substructures[0]

    # instead of raw string equality, canonicalize both SMARTS and compare
    compare_smarts(sub_match.smarts, expected_smarts)
    # atom indices should still be the same
    assert sub_match.idxs_match == expected_idxs_match


def test_detector_stats_count_sites() -> None:
    """Checks that get_stats reports the sites find_sites returned. The
    found, validated, and rejected counters were never incremented."""

    detector = ProtonationSiteDetector()
    _, sites = detector.find_sites(MoleculeRecord("NCCC(=O)O"))
    stats = detector.get_stats()

    assert len(sites) > 0
    assert stats["sites_found"] == len(sites)
    assert stats["sites_validated"] == len(sites)
    assert stats["sites_rejected"] == 0


def test_detector_stats_count_rejected_sites(monkeypatch: pytest.MonkeyPatch) -> None:
    """Checks that sites failing validation are counted as rejected."""

    monkeypatch.setattr(ProtonationSite, "is_valid", lambda self: False)

    detector = ProtonationSiteDetector()
    _, sites = detector.find_sites(MoleculeRecord("NCCC(=O)O"))
    stats = detector.get_stats()

    assert sites == []
    assert stats["sites_found"] > 0
    assert stats["sites_rejected"] == stats["sites_found"]
    assert stats["sites_validated"] == 0
