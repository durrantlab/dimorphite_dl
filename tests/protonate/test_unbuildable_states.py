"""Checks molecules where one site's state has no valid structure."""

import pytest
from rdkit import Chem

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
