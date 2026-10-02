"""Checks molecules where one site's state has no valid structure."""

from typing import NoReturn

import pytest
from rdkit import Chem

from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.protonate.detect import (
    ProtonationSiteDetectionError,
    ProtonationSiteDetector,
)
from dimorphite_dl.protonate.run import Protonate

VERY_ACIDIC_PH = -10000000.0
VERY_BASIC_PH = 10000000.0

# Each has a bridgehead aromatic N that Aromatic_nitrogen_unprotonated
# matches but that cannot carry a +1 charge and still be kekulized. Found by
# running the NCI Diversity Set V.
BRIDGEHEAD_N_MOLECULES = [
    "Cc1ccc2nc(Cl)cc(=O)n2c1",
    "Cn1cnc2ncnn2c1=S",
    "Cc1ccc2c(c1)sc1nncn12",
    "Clc1nnc2c3ccccc3ccn12",
    "O=c1nc(Oc2ccccc2)nc2ncccn12",
]


@pytest.mark.parametrize("smiles", BRIDGEHEAD_N_MOLECULES)
@pytest.mark.parametrize(
    "ph", [VERY_ACIDIC_PH, 7.0, VERY_BASIC_PH], ids=["very_acidic", "7.0", "very_basic"]
)
def test_unbuildable_state_keeps_molecule(smiles: str, ph: float) -> None:
    """Checks that a site with no valid charged structure keeps its parent's
    state instead of producing an invalid variant. The invalid SMILES was
    rejected after enumeration, and when that was the site's only state
    (below about pH 2.3) every variant was rejected and the molecule
    vanished from the output."""

    protonator = Protonate([smiles], ph_min=ph, ph_max=ph, precision=1.0)
    output = protonator.to_list()
    stats = protonator.get_stats()["protonation"]

    assert len(output) > 0, "molecule vanished from the output"
    assert all(Chem.MolFromSmiles(line) is not None for line in output), output
    assert stats["variants_rejected"] == 0, output
    assert stats["fallback_used"] == 0, output


def test_all_variants_rejected_falls_back_to_input(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Checks that a molecule whose every variant fails validation is still
    written once, as its input, rather than dropped. Output is otherwise one
    or more lines per input, and a missing line silently misaligns any
    downstream join."""

    monkeypatch.setattr(Protonate, "_is_smiles_valid", lambda self, smiles: False)

    protonator = Protonate(["CCC(=O)O"])
    assert protonator.to_list() == ["CCC(=O)O"]
    assert protonator.get_stats()["protonation"]["fallback_used"] == 1


def test_detection_error_is_reported_as_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Checks that a crash in site detection is counted as a fallback, not as
    a molecule without sites. The error was swallowed and the molecule was
    reported as having no sites, which hid the failure."""

    def failing_find_sites(
        self: ProtonationSiteDetector, mol_record: MoleculeRecord
    ) -> NoReturn:
        raise ProtonationSiteDetectionError("Detection failed: forced")

    monkeypatch.setattr(ProtonationSiteDetector, "find_sites", failing_find_sites)

    protonator = Protonate(["CCC(=O)O"])
    assert protonator.to_list() == ["CCC(=O)O"]
    stats = protonator.get_stats()["protonation"]
    assert stats["fallback_used"] == 1
    assert stats["molecules_without_sites"] == 0


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        pytest.param("c1nn[nH]n1", "c1nn[n-]n1", id="plain"),
        pytest.param("c1nn[nH:1]n1", "c1nn[n-:1]n1", id="atom_map"),
        pytest.param("c1nn[15nH]n1", "c1nn[15n-]n1", id="isotope"),
    ],
)
def test_labeled_aromatic_nh_can_be_deprotonated(smiles: str, expected: str) -> None:
    """Checks that an aromatic N-H carrying an atom map or isotope label is
    deprotonated like an unlabeled one. The H was removed only when the SMILES
    contained the literal text "[nH-]", so "[nH-:1]" kept its H, could not be
    sanitized, and fell back to the neutral parent at every pH.

    Args:
        smiles: Tetrazole (pKa about 4.9) with or without a label on the N-H.
        expected: The anion, with the label kept.
    """

    output = Protonate([smiles], ph_min=7.4, ph_max=7.4, precision=1.0).to_list()
    assert [Chem.MolToSmiles(Chem.MolFromSmiles(line)) for line in output] == [
        Chem.MolToSmiles(Chem.MolFromSmiles(expected))
    ], output