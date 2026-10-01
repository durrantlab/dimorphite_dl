"""Checks nitrogens that should not fall through to the aliphatic amine site.

The xfail cases are known SMARTS defects in site_substructures.smarts. They are
strict, so fixing a pattern turns the test into an XPASS failure and the
marker has to be removed rather than left in place.
"""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles
from dimorphite_dl.protonate.run import Protonate

VERY_ACIDIC_PH = -10000000.0


def canonical_set(smiles_list: list[str]) -> set[str]:
    """Canonicalize outputs so the comparison ignores SMILES writing order.

    Args:
        smiles_list: SMILES strings returned by protonate_smiles.

    Returns:
        Set of canonical isomeric SMILES.
    """
    return {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smiles_list}


@pytest.mark.parametrize(
    "smiles",
    [
        # Control: CH2 next to N already matches Anilines_secondary.
        pytest.param("CCNc1ccccc1", id="N-ethylaniline"),
        pytest.param("CC(C)Nc1ccccc1", id="N-isopropylaniline"),
        pytest.param("CC(C)N(C)c1ccccc1", id="N-isopropyl-N-methylaniline"),
        pytest.param("CN(C)S(C)(=O)=O", id="tertiary_sulfonamide"),
        pytest.param("CN(C)C(=O)OC(C)(C)C", id="boc_carbamate"),
    ],
)
def test_nitrogen_stays_neutral_at_physiological_ph(smiles: str) -> None:
    # A detection or protonation error falls back to echoing the input, which
    # would also satisfy the first assertion, so fallbacks are checked too.
    protonator = Protonate([smiles], ph_min=7.4, ph_max=7.4)
    output = protonator.to_list()
    assert canonical_set(output) == canonical_set([smiles]), output
    assert protonator.get_stats()["protonation"]["fallback_used"] == 0


@pytest.mark.parametrize(
    ("smiles", "protonated"),
    [
        pytest.param("CCNc1ccccc1", "CC[NH2+]c1ccccc1", id="N-ethylaniline"),
        pytest.param("CC(C)Nc1ccccc1", "CC(C)[NH2+]c1ccccc1", id="N-isopropylaniline"),
        pytest.param(
            "CC(C)N(C)c1ccccc1",
            "CC(C)[NH+](C)c1ccccc1",
            id="N-isopropyl-N-methylaniline",
        ),
    ],
)
def test_aniline_nitrogen_is_still_a_site(smiles: str, protonated: str) -> None:
    # Neutral at pH 7.4 is also what a molecule with no site at all would give.
    # Protonation at very low pH shows the nitrogen is still detected, and the
    # neutral result at 7.4 shows it was assigned an aniline pKa, not an amine pKa.
    output = protonate_smiles(
        smiles, ph_min=VERY_ACIDIC_PH, ph_max=VERY_ACIDIC_PH, precision=0.5
    )
    assert canonical_set(output) == canonical_set([protonated]), output
