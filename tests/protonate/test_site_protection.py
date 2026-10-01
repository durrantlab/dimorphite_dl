"""Checks that site detection gives each ionizable group exactly its own
sites, without one group hiding or duplicating another's."""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

VERY_ACIDIC_PH = -10000000.0
VERY_BASIC_PH = 10000000.0


def canonical(smiles: str) -> str:
    """Canonicalizes so expected values can be written in any valid form.

    Args:
        smiles: A valid SMILES string.

    Returns:
        RDKit's canonical isomeric SMILES for smiles.
    """

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


def protonate_labeled(smiles: str, ph: float) -> list[tuple[str, list[str]]]:
    """Runs the pipeline at a single pH with state labels, since a duplicated
    or missing site shows up in the labels even when the SMILES looks right.

    Args:
        smiles: The input SMILES string.
        ph: Used as both the minimum and maximum pH.

    Returns:
        One (canonical SMILES, labels) pair per output state.
    """

    output = protonate_smiles(
        smiles, ph_min=ph, ph_max=ph, precision=0.5, label_states=True
    )
    result = []
    for line in output:
        smiles_out, _, states = line.partition(",")
        result.append((canonical(smiles_out), states.split("\t") if states else []))
    return result


SHARED_CONTEXT_MOLECULES = [
    # [input, very acidic, very basic, number of sites, id]
    [
        "NCP(=O)(O)O",
        "[NH3+]CP(=O)(O)O",
        "NCP(=O)([O-])[O-]",
        3,
        "phosphonate_carbon_shared_with_amine",
    ],
    [
        "O=C(O)c1ccc[nH]1",
        "O=C(O)c1ccc[nH]1",
        "O=C([O-])c1ccc[n-]1",
        2,
        "carboxyl_carbon_shared_with_pyrrole",
    ],
    ["OP(=O)(O)O", "O=P(O)(O)O", "O=P([O-])([O-])O", 2, "phosphoric_acid"],
    [
        "OP(=O)(O)OP(=O)(O)O",
        "O=P(O)(O)OP(=O)(O)O",
        "O=P([O-])([O-])OP(=O)([O-])[O-]",
        4,
        "pyrophosphate",
    ],
]


@pytest.mark.parametrize(
    ("smiles", "acidic", "basic", "n_sites"),
    [group[:4] for group in SHARED_CONTEXT_MOLECULES],
    ids=[group[4] for group in SHARED_CONTEXT_MOLECULES],
)
@pytest.mark.parametrize(
    ("ph", "state"),
    [(VERY_ACIDIC_PH, "PROTONATED"), (VERY_BASIC_PH, "DEPROTONATED")],
    ids=["very_acidic", "very_basic"],
)
def test_shared_context_sites(
    smiles: str, acidic: str, basic: str, n_sites: int, ph: float, state: str
) -> None:
    """Checks that a carbon one group used only as its attachment point does
    not hide a later group's site, and that one group's repeated matches do
    not assign extra sites. Every matched atom was locked, carbons included,
    and all matches of a pattern were accepted before any was locked, so the
    amine of NCP(=O)(O)O was never protonated and phosphoric acid lost all
    three protons."""

    expected = acidic if state == "PROTONATED" else basic
    assert protonate_labeled(smiles, ph) == [(canonical(expected), [state] * n_sites)]


def test_primary_amide_is_one_site() -> None:
    """Checks that the two N-H matches of a primary amide give one site. Both
    were accepted, so the same nitrogen was listed, and charged, twice."""

    assert protonate_labeled("CC(N)=O", VERY_BASIC_PH) == [
        (canonical("CC([NH-])=O"), ["DEPROTONATED"])
    ]
