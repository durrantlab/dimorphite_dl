import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

VERY_BASIC_PH = 10000000.0


def canonical_set(smiles_list: list[str]) -> set[str]:
    """Compare outputs by structure rather than by SMILES spelling.

    Args:
        smiles_list: SMILES strings to canonicalize.

    Returns:
        The set of canonical SMILES.
    """
    return {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smiles_list}


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("Cc1c[nH]cn1", id="imidazole"),
        pytest.param("c1cn[nH]c1", id="pyrazole"),
        pytest.param("c1nc[nH]n1", id="1,2,4-triazole"),
        pytest.param("c1c[nH]nn1", id="1,2,3-triazole"),
        pytest.param("c1ccc2[nH]cnc2c1", id="benzimidazole"),
        pytest.param("Nc1ncnc2[nH]cnc12", id="adenine"),
    ],
)
def test_azole_nh_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Azole N-H pKas are far above 7.4, so no anion should be produced.

    Args:
        smiles: A neutral azole.
    """
    output = protonate_smiles(smiles, ph_min=7.4, ph_max=7.4)
    assert canonical_set(output) == canonical_set([smiles]), output


def test_imidazole_nh_deprotonates_at_very_basic_ph() -> None:
    """Confirm the Azole_NH site is still wired to the ring nitrogen."""
    output = protonate_smiles(
        "c1c[nH]cn1", ph_min=VERY_BASIC_PH, ph_max=VERY_BASIC_PH, precision=0.5
    )
    assert canonical_set(output) == canonical_set(["c1c[n-]cn1"]), output