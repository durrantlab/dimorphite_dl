"""Checks the site entries whose pKas come from the literature rather than the
training data: Hydrazoic_acid, Carboxamide, and Diazine."""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

PHYSIOLOGICAL_PH = 7.4

# Between the Diazine window (about 0.7 to 2.8) and the upper end of the
# Aromatic_nitrogen_unprotonated window (about 2.3 to 6.4), so the two rules
# give different answers here.
DIAZINE_CONTROL_PH = 4.0


def canonical_set(smiles_list: list[str]) -> set[str]:
    """Compare outputs by structure rather than by SMILES spelling.

    Args:
        smiles_list: SMILES strings to canonicalize.

    Returns:
        The set of canonical SMILES.
    """
    return {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smiles_list}


def protonate_at(smiles: str, ph: float) -> list[str]:
    """Runs a single pH so each test reads as one expected state.

    Args:
        smiles: The input SMILES string.
        ph: The pH to protonate at.

    Returns:
        The output SMILES.
    """
    return protonate_smiles(smiles, ph_min=ph, ph_max=ph)


def test_hydrazoic_acid_is_azide_at_physiological_ph() -> None:
    """HN3 has a pKa of 4.65. Azide alone stopped at neutral HN3, because the
    neutralizer adds a second proton that Azide cannot remove."""
    output = protonate_at("[N-]=[N+]=N", PHYSIOLOGICAL_PH)
    assert canonical_set(output) == canonical_set(["[N-]=[N+]=[N-]"]), output


def test_organic_azide_is_unchanged() -> None:
    """Hydrazoic_acid needs an H on both terminal nitrogens, so R-N3 must
    still go to Azide."""
    output = protonate_at("CN=[N+]=[N-]", PHYSIOLOGICAL_PH)
    assert canonical_set(output) == canonical_set(["CN=[N+]=[N-]"]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(N)=O", id="acetamide"),
        pytest.param("CNC(C)=O", id="N-methylacetamide"),
        pytest.param("CC(=O)Nc1ccccc1", id="acetanilide"),
        pytest.param("NC(=O)c1ccccc1", id="benzamide"),
        pytest.param("O=C1CCCCCN1", id="caprolactam"),
        pytest.param("CCN(CC)CC(=O)Nc1c(C)cccc1C", id="lidocaine"),
        pytest.param("NCC(=O)NCC(=O)O", id="glycylglycine"),
    ],
)
def test_carboxamide_nh_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Simple amide N-H has a pKa near 15, so no amide anion should appear.
    Other sites (e.g., an amine) may still vary, so only N- is checked.

    Args:
        smiles: A molecule with a simple carboxamide N-H.
    """
    for out in protonate_at(smiles, PHYSIOLOGICAL_PH):
        mol = Chem.MolFromSmiles(out)
        charged_n = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "N" and atom.GetFormalCharge() < 0
        ]
        assert not charged_n, out


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("c1cncnc1", id="pyrimidine"),
        pytest.param("c1cnccn1", id="pyrazine"),
        pytest.param("c1ccnnc1", id="pyridazine"),
        pytest.param("c1ccc2ncncc2c1", id="quinazoline"),
        pytest.param("c1ccc2nccnc2c1", id="quinoxaline"),
        pytest.param("c1ncncn1", id="1,3,5-triazine"),
        pytest.param("NC(=O)c1cnccn1", id="pyrazinamide"),
    ],
)
def test_diazine_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Diazine pKaHs are 0.6 to 3.5, so neither a cation nor a dication
    should appear.

    Args:
        smiles: A neutral diazine.
    """
    output = protonate_at(smiles, PHYSIOLOGICAL_PH)
    assert canonical_set(output) == canonical_set([smiles]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("Nc1ncccn1", id="2-aminopyrimidine"),
        pytest.param("COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", id="trimethoprim"),
    ],
)
def test_donor_substituted_diazine_keeps_generic_rule(smiles: str) -> None:
    """Amino groups raise diazine pKaH sharply, so these rings are excluded
    from Diazine. A protonated variant at pH 4 shows they still get the
    Aromatic_nitrogen_unprotonated pKa.

    Args:
        smiles: An amino-substituted pyrimidine.
    """
    output = protonate_at(smiles, DIAZINE_CONTROL_PH)
    assert any("+" in out for out in output), output
