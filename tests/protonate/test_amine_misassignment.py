"""Checks nitrogens that should not fall through to the aliphatic amine site.

The xfail cases are known SMARTS defects in site_substructures.smarts. They are
strict, so fixing a pattern turns the test into an XPASS failure and the
marker has to be removed rather than left in place.
"""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

ANILINE_H_COUNT = pytest.mark.xfail(
    strict=True,
    reason="Anilines_secondary/tertiary use [!H] (H-count query) for 'heavy "
    "atom', so an N-CH(R)R' substituent falls through to the amine site",
)
ACYL_TERTIARY_N = pytest.mark.xfail(
    strict=True,
    reason="Amines_primary_secondary_tertiary matches N-acyl and N-sulfonyl "
    "nitrogens that lack an N-H, so the Amide and Sulfonamide patterns miss them",
)


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
        pytest.param("CC(C)Nc1ccccc1", id="N-isopropylaniline", marks=ANILINE_H_COUNT),
        pytest.param(
            "CC(C)N(C)c1ccccc1",
            id="N-isopropyl-N-methylaniline",
            marks=ANILINE_H_COUNT,
        ),
        pytest.param(
            "CN(C)S(C)(=O)=O", id="tertiary_sulfonamide", marks=ACYL_TERTIARY_N
        ),
        pytest.param("CN(C)C(=O)OC(C)(C)C", id="boc_carbamate", marks=ACYL_TERTIARY_N),
    ],
)
def test_nitrogen_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Checks that a weakly basic nitrogen gives only the neutral form at pH
    7.4. Misassigning it to the aliphatic amine site (pKa near 8) adds a
    spurious cation and doubles the variant count."""
    output = protonate_smiles(smiles, ph_min=7.4, ph_max=7.4)

    assert canonical_set(output) == canonical_set([smiles])
