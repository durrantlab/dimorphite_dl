"""Checks that the output does not depend on how charged the input is."""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

# Wide enough that nearly every site is BOTH at pH 7, so every state of each
# group is compared rather than only the dominant one.
WIDE_PRECISION = 10.0

# [neutral input, protonated form, deprotonated form, group]. The charged
# forms are the expected states from test_pka_average in test_run.py.
GROUPS = [
    ["C#CCO", "C#CCO", "C#CC[O-]", "Alcohol"],
    ["C(=O)N", "NC=O", "[NH-]C=O", "Amide"],
    ["CC(=O)NOC(C)=O", "CC(=O)NOC(C)=O", "CC(=O)[N-]OC(C)=O", "Amide_electronegative"],
    ["COC(=N)N", "COC(N)=[NH2+]", "COC(=N)N", "AmidineGuanidine2"],
    [
        "Brc1ccc(C2NCCS2)cc1",
        "Brc1ccc(C2[NH2+]CCS2)cc1",
        "Brc1ccc(C2NCCS2)cc1",
        "Amines_primary_secondary_tertiary",
    ],
    [
        "CC(=O)[n+]1ccc(N)cc1",
        "CC(=O)[n+]1ccc([NH3+])cc1",
        "CC(=O)[n+]1ccc(N)cc1",
        "Anilines_primary",
    ],
    ["CCNc1ccccc1", "CC[NH2+]c1ccccc1", "CCNc1ccccc1", "Anilines_secondary"],
    ["Cc1ccccc1N(C)C", "Cc1ccccc1[NH+](C)C", "Cc1ccccc1N(C)C", "Anilines_tertiary"],
    [
        "BrC1=CC2=C(C=C1)NC=C2",
        "Brc1ccc2[nH]ccc2c1",
        "Brc1ccc2[n-]ccc2c1",
        "Indole_pyrrole",
    ],
    ["O=c1cc[nH]cc1", "O=c1cc[nH]cc1", "O=c1cc[n-]cc1", "Aromatic_nitrogen_protonated"],
    ["C-N=[N+]=[N@H]", "CN=[N+]=N", "CN=[N+]=[N-]", "Azide"],
    ["BrC(C(O)=O)CBr", "O=C(O)C(Br)CBr", "O=C([O-])C(Br)CBr", "Carboxyl"],
    ["NC(NN=O)=N", "NC(=[NH2+])NN=O", "N=C(N)NN=O", "AmidineGuanidine1"],
    [
        "C(F)(F)(F)C(=O)NC(=O)C",
        "CC(=O)NC(=O)C(F)(F)F",
        "CC(=O)[N-]C(=O)C(F)(F)F",
        "Imide",
    ],
    ["O=C(C)NC(C)=O", "CC(=O)NC(C)=O", "CC(=O)[N-]C(C)=O", "Imide2"],
    [
        "CC(C)(C)C(N(C)O)=O",
        "CN(O)C(=O)C(C)(C)C",
        "CN([O-])C(=O)C(C)(C)C",
        "N-hydroxyamide",
    ],
    ["C[N+](O)=O", "C[N+](=O)O", "C[N+](=O)[O-]", "Nitro"],
    ["O=C1C=C(O)CC1", "O=C1C=C(O)CC1", "O=C1C=C([O-])CC1", "O=C-C=C-OH"],
    ["C1CC1OO", "OOC1CC1", "[O-]OC1CC1", "Peroxide2"],
    ["C(=O)OO", "O=COO", "O=CO[O-]", "Peroxide1"],
    ["Brc1cc(O)cc(Br)c1", "Oc1cc(Br)cc(Br)c1", "[O-]c1cc(Br)cc(Br)c1", "Phenol"],
    ["CC(=O)c1ccc(S)cc1", "CC(=O)c1ccc(S)cc1", "CC(=O)c1ccc([S-])cc1", "Phenyl_Thiol"],
    [
        "C=CCOc1ccc(C(=O)O)cc1",
        "C=CCOc1ccc(C(=O)O)cc1",
        "C=CCOc1ccc(C(=O)[O-])cc1",
        "Phenyl_carboxyl",
    ],
    ["COP(=O)(O)OC", "COP(=O)(O)OC", "COP(=O)([O-])OC", "Phosphate_diester"],
    ["CP(C)(=O)O", "CP(C)(=O)O", "CP(C)(=O)[O-]", "Phosphinic_acid"],
    [
        "CC(C)OP(C)(=O)O",
        "CC(C)OP(C)(=O)O",
        "CC(C)OP(C)(=O)[O-]",
        "Phosphonate_ester",
    ],
    ["CC1(C)OC(=O)NC1=O", "CC1(C)OC(=O)NC1=O", "CC1(C)OC(=O)[N-]C1=O", "Ringed_imide1"],
    ["O=C(N1)C=CC1=O", "O=C1C=CC(=O)N1", "O=C1C=CC(=O)[N-]1", "Ringed_imide2"],
    ["O=S(OC)(O)=O", "COS(=O)(=O)O", "COS(=O)(=O)[O-]", "Sulfate"],
    [
        "COc1ccc(S(=O)O)cc1",
        "COc1ccc(S(=O)O)cc1",
        "COc1ccc(S(=O)[O-])cc1",
        "Sulfinic_acid",
    ],
    ["CS(N)(=O)=O", "CS(N)(=O)=O", "CS([NH-])(=O)=O", "Sulfonamide"],
    [
        "CC(=O)CSCCS(O)(=O)=O",
        "CC(=O)CSCCS(=O)(=O)O",
        "CC(=O)CSCCS(=O)(=O)[O-]",
        "Sulfonate",
    ],
    ["CC(=O)S", "CC(=O)S", "CC(=O)[S-]", "Thioic_acid"],
    ["C(C)(C)(C)(S)", "CC(C)(C)S", "CC(C)(C)[S-]", "Thiol"],
    [
        "Brc1cc[nH+]cc1",
        "Brc1cc[nH+]cc1",
        "Brc1ccncc1",
        "Aromatic_nitrogen_unprotonated",
    ],
    [
        "C=C(O)c1c(C)cc(C)cc1C",
        "C=C(O)c1c(C)cc(C)cc1C",
        "C=C([O-])c1c(C)cc(C)cc1C",
        "Vinyl_alcohol",
    ],
    ["CC(=O)ON", "CC(=O)O[NH3+]", "CC(=O)ON", "Primary_hydroxyl_amine"],
    [
        "O=P(O)(O)OCCCC",
        "CCCCOP(=O)([O-])O",
        "CCCCOP(=O)([O-])[O-]",
        "Phosphate",
    ],
    ["CC(P(O)(O)=O)C", "CC(C)P(=O)([O-])O", "CC(C)P(=O)([O-])[O-]", "Phosphonate"],
]


def protonate_normalized(smiles: str) -> list[str]:
    """Runs the pipeline at pH 7 with labels and canonicalizes each output
    SMILES, so forms written differently compare equal.

    Args:
        smiles: The input SMILES string.

    Returns:
        Sorted output lines, each a canonical SMILES plus its labels.
    """

    output = protonate_smiles(
        smiles, ph_min=7.0, ph_max=7.0, precision=WIDE_PRECISION, label_states=True
    )
    lines = []
    for line in output:
        smiles_out, sep, states = line.partition(",")
        mol = Chem.MolFromSmiles(smiles_out)
        lines.append(Chem.MolToSmiles(mol, isomericSmiles=True) + sep + states)
    return sorted(lines)


@pytest.mark.parametrize(
    ("neutral", "protonated", "deprotonated"),
    [group[:3] for group in GROUPS],
    ids=[group[3] for group in GROUPS],
)
def test_output_does_not_depend_on_input_charge_state(
    neutral: str, protonated: str, deprotonated: str
) -> None:
    """Checks that each group gives the same states and labels whether it is
    written neutral, protonated, or deprotonated. Neutralization handles
    only the charged forms it has a rule for, so a missing rule shows up here
    as a charged input that keeps its charge or loses its site."""

    expected = protonate_normalized(neutral)

    for charged_form in [protonated, deprotonated]:
        assert protonate_normalized(charged_form) == expected, charged_form
