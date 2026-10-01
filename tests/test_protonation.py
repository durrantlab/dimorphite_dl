"""Checks the protonation states Dimorphite-DL assigns to each ionizable
group, alone and in combination."""

import os
import sys
from io import StringIO
from typing import Callable, Dict, List, Tuple, Union

import pytest
from rdkit import Chem

import dimorphite_dl

SMARTS_FILE = os.path.join(
    os.path.dirname(os.path.realpath(dimorphite_dl.__file__)),
    "site_substructures.smarts",
)

VERY_ACIDIC_PH = -10000000.0
VERY_BASIC_PH = 10000000.0
DEFAULT_PKA_PRECISION = 0.5

SINGLE_SITE_GROUPS = [
    # [input smiles, protonated, deprotonated, category]
    ["C#CCO", "C#CCO", "C#CC[O-]", "Alcohol"],
    ["C(=O)N", "NC=O", "[NH-]C=O", "Amide"],
    [
        "CC(=O)NOC(C)=O",
        "CC(=O)NOC(C)=O",
        "CC(=O)[N-]OC(C)=O",
        "Amide_electronegative",
    ],
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
    [
        "Cc1ccccc1N(C)C",
        "Cc1ccccc1[NH+](C)C",
        "Cc1ccccc1N(C)C",
        "Anilines_tertiary",
    ],
    [
        "BrC1=CC2=C(C=C1)NC=C2",
        "Brc1ccc2[nH]ccc2c1",
        "Brc1ccc2[n-]ccc2c1",
        "Indole_pyrrole",
    ],
    # N-protonated 4-pyridone has no valid structure (the real cation is
    # O-protonated), so the protonated state is the neutral pyridone.
    [
        "BrC1=CNC=C(C1=O)Br",
        "O=c1c(Br)c[nH]cc1Br",
        "O=c1c(Br)c[nH]cc1Br",
        "Aromatic_nitrogen_protonated",
    ],
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
    [
        "Brc1cc(O)cc(Br)c1",
        "Oc1cc(Br)cc(Br)c1",
        "[O-]c1cc(Br)cc(Br)c1",
        "Phenol",
    ],
    [
        "CC(=O)c1ccc(S)cc1",
        "CC(=O)c1ccc(S)cc1",
        "CC(=O)c1ccc([S-])cc1",
        "Phenyl_Thiol",
    ],
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
    [
        "CC1(C)OC(=O)NC1=O",
        "CC1(C)OC(=O)NC1=O",
        "CC1(C)OC(=O)[N-]C1=O",
        "Ringed_imide1",
    ],
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
]


PHOSPHORUS_GROUPS = [
    # [input smiles, protonated, singly deprotonated, deprotonated, category]
    [
        "O=P(O)(O)OCCCC",
        "CCCCOP(=O)(O)O",
        "CCCCOP(=O)([O-])O",
        "CCCCOP(=O)([O-])[O-]",
        "Phosphate",
    ],
    [
        "CC(P(O)(O)=O)C",
        "CC(C)P(=O)(O)O",
        "CC(C)P(=O)([O-])O",
        "CC(C)P(=O)([O-])[O-]",
        "Phosphonate",
    ],
]

# [input smiles, protonated, deprotonated]. Written in any valid form; the
# tests canonicalize them.
MULTIPLE_SITE_MOLECULES = [
    ["NCCCC(=O)O", "[NH3+]CCCC(=O)O", "NCCCC(=O)[O-]"],
    ["NCc1ccc(O)cc1", "[NH3+]Cc1ccc(O)cc1", "NCc1ccc([O-])cc1"],
    [
        "OC(=O)CCC(N)C(=O)O",
        "OC(=O)CCC([NH3+])C(=O)O",
        "[O-]C(=O)CCC(N)C(=O)[O-]",
    ],
    ["NCc1ccc2[nH]ccc2c1", "[NH3+]Cc1ccc2[nH]ccc2c1", "NCc1ccc2[n-]ccc2c1"],
]

# [input smiles, protonated, deprotonated]. Each has two sites whose matches
# share context atoms, or one pattern with several matches on the same sites.
SHARED_CONTEXT_MOLECULES = [
    # The phosphonate's carbon is also the amine's only carbon.
    ["NCP(=O)(O)O", "[NH3+]CP(=O)(O)O", "NCP(=O)([O-])[O-]"],
    # The carboxyl's ring carbon is also part of the pyrrole match.
    ["O=C(O)c1ccc[nH]1", "O=C(O)c1ccc[nH]1", "O=C([O-])c1ccc[n-]1"],
    # Three interchangeable OH groups; only two are sites.
    ["OP(=O)(O)O", "O=P(O)(O)O", "O=P([O-])([O-])O"],
    # Both phosphates use the same bridging O as context.
    [
        "OP(=O)(O)OP(=O)(O)O",
        "O=P(O)(O)OP(=O)(O)O",
        "O=P([O-])([O-])OP(=O)([O-])[O-]",
    ],
]

# A precision this wide makes every site BOTH at any pH.
EVERY_SITE_BOTH_PRECISION = 1000000.0

# Eight carboxyl groups on distinct, non-equivalent carbons, so every one of
# the 2**8 charge patterns is a different molecule.
EIGHT_CARBOXYLS = "CC(C(=O)O)" * 8
EIGHT_CARBOXYL_STATES = 2**8


def group_id(group: List[str]) -> str:
    """Names each parametrized case after its category so failures are easy
    to find in pytest's output.

    Args:
        group: A row of SINGLE_SITE_GROUPS or PHOSPHORUS_GROUPS.

    Returns:
        The category name, the last entry of the row.
    """

    return group[-1]


@pytest.fixture(scope="module")
def average_pkas() -> Dict[str, List[float]]:
    """Reads each group's mean pKa values from the shipped SMARTS file, so the
    at-pKa tests follow the parameters rather than hard-coding them.

    Returns:
        Group name (with any "*" removed) mapped to its mean pKa per site.
    """

    pkas: Dict[str, List[float]] = {}
    with open(SMARTS_FILE) as smarts:
        for line in smarts:
            splits = line.split()
            if len(splits) == 0:
                continue

            # Columns: name, SMARTS, then (site, mean, std) triples.
            pkas[splits[0].replace("*", "")] = [float(x) for x in splits[3::3]]
    return pkas


def protonate(
    smiles: str,
    ph: float,
    pka_precision: float,
    max_variants: int = dimorphite_dl.DEFAULT_MAX_VARIANTS,
) -> List[List[str]]:
    """Runs Dimorphite-DL on one molecule at a single pH, with state labels,
    the way the command line would.

    Args:
        smiles: The input SMILES string.
        ph: Used as both the minimum and maximum pH.
        pka_precision: Number of standard deviations to consider.
        max_variants: The most states to return.

    Returns:
        One [smiles, first label, other labels...] list per output state.
    """

    args: Dict[str, Union[str, float, bool]] = {
        "min_ph": ph,
        "max_ph": ph,
        "pka_precision": pka_precision,
        "max_variants": max_variants,
        "smiles": smiles,
        "label_states": True,
    }
    return [line.split() for line in dimorphite_dl.Protonate(args)]


def normalize_smiles(smiles: str) -> str:
    """Canonicalizes with the installed RDKit, so comparisons do not depend on
    the RDKit version that wrote the expected strings. Unparseable strings
    cannot be canonicalized and are compared as is.

    Args:
        smiles: A SMILES string, valid or not.

    Returns:
        The canonical isomeric SMILES, or smiles itself if it does not parse.
    """

    mol = Chem.MolFromSmiles(smiles)
    return smiles if mol is None else Chem.MolToSmiles(mol, isomericSmiles=True)


def check_protonation(
    smiles: str,
    ph: float,
    expected_smiles: List[str],
    labels: List[str],
    pka_precision: float = DEFAULT_PKA_PRECISION,
) -> List[List[str]]:
    """Asserts that a molecule yields exactly the expected states at a pH.

    Every output must also be a parseable SMILES string, so that bad hydrogen
    bookkeeping (e.g., [nH-]) is caught.

    Args:
        smiles: The input SMILES string.
        ph: Used as both the minimum and maximum pH.
        expected_smiles: The canonical SMILES of every expected state.
            Repeats count once, since Dimorphite-DL drops duplicate states.
        labels: The allowed labels. Every site's label is checked, not just
            the first.
        pka_precision: Number of standard deviations to consider.

    Returns:
        The output, so callers can make further assertions (e.g., on the
        order of the site labels).
    """

    output = protonate(smiles, ph, pka_precision)
    output_smiles = [line[0] for line in output]

    invalid = [s for s in output_smiles if Chem.MolFromSmiles(s) is None]
    assert invalid == [], "invalid SMILES produced: " + str(invalid)

    # Lists rather than sets, so a repeated output cannot stand in for a
    # missing expected state.
    expected = sorted(set(normalize_smiles(s) for s in expected_smiles))
    assert sorted(normalize_smiles(s) for s in output_smiles) == expected, output
    assert set(label for line in output for label in line[1:]) <= set(labels), output
    return output


def test_check_protonation_rejects_repeated_state_for_missing_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Checks that the helper fails when one expected state is replaced by a
    second copy of another. It compared only the count and a subset, so
    [A, A] passed for an expectation of {A, B}."""

    monkeypatch.setattr(
        sys.modules[__name__],
        "protonate",
        lambda *args: [["CCO", "BOTH"], ["CCO", "BOTH"]],
    )

    with pytest.raises(AssertionError):
        check_protonation("CCO", 7.0, ["CCO", "CC[O-]"], ["BOTH"])


@pytest.mark.parametrize("group", SINGLE_SITE_GROUPS, ids=group_id)
def test_single_site_very_acidic(group: List[str]) -> None:
    """Checks that every group is fully protonated at extremely low pH."""

    smiles, protonated, deprotonated, category = group
    check_protonation(smiles, VERY_ACIDIC_PH, [protonated], ["PROTONATED"])


@pytest.mark.parametrize("group", SINGLE_SITE_GROUPS, ids=group_id)
def test_single_site_very_basic(group: List[str]) -> None:
    """Checks that every group is fully deprotonated at extremely high pH."""

    smiles, protonated, deprotonated, category = group
    check_protonation(smiles, VERY_BASIC_PH, [deprotonated], ["DEPROTONATED"])


@pytest.mark.parametrize("group", SINGLE_SITE_GROUPS, ids=group_id)
def test_single_site_at_category_pka(
    group: List[str], average_pkas: Dict[str, List[float]]
) -> None:
    """Checks that both states are produced at a group's own mean pKa."""

    smiles, protonated, deprotonated, category = group
    check_protonation(
        smiles, average_pkas[category][0], [protonated, deprotonated], ["BOTH"]
    )


@pytest.mark.parametrize("group", PHOSPHORUS_GROUPS, ids=group_id)
def test_phosphorus_very_acidic(group: List[str]) -> None:
    """Checks that both acidic sites stay protonated at extremely low pH."""

    smiles, protonated, mix, deprotonated, category = group
    check_protonation(smiles, VERY_ACIDIC_PH, [protonated], ["PROTONATED"])


@pytest.mark.parametrize("group", PHOSPHORUS_GROUPS, ids=group_id)
def test_phosphorus_very_basic(group: List[str]) -> None:
    """Checks that both acidic sites are deprotonated at extremely high pH."""

    smiles, protonated, mix, deprotonated, category = group
    check_protonation(smiles, VERY_BASIC_PH, [deprotonated], ["DEPROTONATED"])


@pytest.mark.parametrize("group", PHOSPHORUS_GROUPS, ids=group_id)
def test_phosphorus_at_first_pka(
    group: List[str], average_pkas: Dict[str, List[float]]
) -> None:
    """Checks that only the first site is ambiguous at the first pKa."""

    smiles, protonated, mix, deprotonated, category = group
    output = check_protonation(
        smiles,
        average_pkas[category][0],
        [mix, protonated],
        ["BOTH", "PROTONATED"],
    )

    # Sites are listed in SMARTS order: the pKa1 site, then the pKa2 site.
    assert all(line[1:] == ["BOTH", "PROTONATED"] for line in output), output


@pytest.mark.parametrize("group", PHOSPHORUS_GROUPS, ids=group_id)
def test_phosphorus_at_second_pka(
    group: List[str], average_pkas: Dict[str, List[float]]
) -> None:
    """Checks that the first site is deprotonated and only the second is
    ambiguous at the second pKa."""

    smiles, protonated, mix, deprotonated, category = group
    output = check_protonation(
        smiles,
        average_pkas[category][1],
        [mix, deprotonated],
        ["DEPROTONATED", "BOTH"],
    )

    # Sites are listed in SMARTS order: the pKa1 site, then the pKa2 site.
    assert all(line[1:] == ["DEPROTONATED", "BOTH"] for line in output), output


@pytest.mark.parametrize("group", PHOSPHORUS_GROUPS, ids=group_id)
def test_phosphorus_between_pkas(
    group: List[str], average_pkas: Dict[str, List[float]]
) -> None:
    """Checks that a wide precision between the two pKas yields all three
    states, with the duplicate singly deprotonated state removed."""

    smiles, protonated, mix, deprotonated, category = group
    ph = 0.5 * (average_pkas[category][0] + average_pkas[category][1])
    check_protonation(
        smiles, ph, [mix, deprotonated, protonated], ["BOTH"], pka_precision=5.0
    )


def test_check_protonation_verifies_every_site_label() -> None:
    """Checks that check_protonation fails when only a later site's label is
    wrong. At pH 6 the carboxyl (listed first) is DEPROTONATED and the amine
    is PROTONATED, so allowing only the first label must fail."""

    zwitterion = "[NH3+]CCCC(=O)[O-]"
    with pytest.raises(AssertionError):
        check_protonation("NCCCC(=O)O", 6.0, [zwitterion], ["DEPROTONATED"])

    check_protonation("NCCCC(=O)O", 6.0, [zwitterion], ["DEPROTONATED", "PROTONATED"])


@pytest.mark.parametrize(
    "molecule", MULTIPLE_SITE_MOLECULES, ids=lambda molecule: molecule[0]
)
def test_multiple_sites_very_acidic(
    molecule: List[str], canonical_smiles: Callable[[str], str]
) -> None:
    """Checks molecules whose first charged site can reorder the canonical
    atoms that later site indices refer to."""

    smiles, protonated, deprotonated = molecule
    check_protonation(
        smiles, VERY_ACIDIC_PH, [canonical_smiles(protonated)], ["PROTONATED"]
    )


@pytest.mark.parametrize(
    "molecule", MULTIPLE_SITE_MOLECULES, ids=lambda molecule: molecule[0]
)
def test_multiple_sites_very_basic(
    molecule: List[str], canonical_smiles: Callable[[str], str]
) -> None:
    """Same as the acidic case. The indole also checks that deprotonating a
    bracketed [nH] does not make the molecule disappear from the output."""

    smiles, protonated, deprotonated = molecule
    check_protonation(
        smiles, VERY_BASIC_PH, [canonical_smiles(deprotonated)], ["DEPROTONATED"]
    )


@pytest.mark.parametrize(
    "molecule", SHARED_CONTEXT_MOLECULES, ids=lambda molecule: molecule[0]
)
@pytest.mark.parametrize(
    "ph, state",
    [(VERY_ACIDIC_PH, "PROTONATED"), (VERY_BASIC_PH, "DEPROTONATED")],
    ids=["very_acidic", "very_basic"],
)
def test_shared_context_sites(
    molecule: List[str],
    ph: float,
    state: str,
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that an atom used only as context by one group does not hide a
    later group's site, and that one group's repeated matches do not assign
    extra sites."""

    smiles, protonated, deprotonated = molecule
    expected = protonated if state == "PROTONATED" else deprotonated
    check_protonation(smiles, ph, [canonical_smiles(expected)], [state])


def test_bracketed_atom_gains_hydrogen_on_protonation() -> None:
    """Checks that protonating an atom written in brackets adds its H.
    Bracketed atoms get no implicit Hs, so without this the isotope-labeled
    amine came out as the invalid [15NH2+]."""

    check_protonation("[15NH2]CC", VERY_ACIDIC_PH, ["CC[15NH3+]"], ["PROTONATED"])


def test_invalid_site_state_keeps_other_sites(
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that a site with no valid protonated structure stays neutral
    without undoing the protonation of the molecule's other sites."""

    check_protonation(
        "O=c1c(Br)c[nH]cc1CN",
        VERY_ACIDIC_PH,
        [canonical_smiles("O=c1c(Br)c[nH]cc1C[NH3+]")],
        ["PROTONATED"],
    )


@pytest.mark.parametrize(
    "smiles",
    ["C[N+](C)(C)[O-]", "[O-][n+]1ccccc1", "CC=[N+](C)[O-]"],
    ids=["amine_oxide", "pyridine_n_oxide", "nitrone"],
)
def test_n_oxide_keeps_oxide(smiles: str) -> None:
    """Checks that neutralization does not turn an N-oxide into an N-hydroxy
    cation. No site restores the O-, so the trimethylamine oxide came out as
    [N+]-OH, and the pyridine N-oxide was misread as a phenol."""

    check_protonation(smiles, 7.0, [smiles], [])


@pytest.mark.parametrize(
    "ph, expected, state",
    [
        (VERY_ACIDIC_PH, "C[N+](=O)O", "PROTONATED"),
        (VERY_BASIC_PH, "C[N+](=O)[O-]", "DEPROTONATED"),
    ],
    ids=["very_acidic", "very_basic"],
)
def test_charged_nitro_input_still_handled(
    ph: float, expected: str, state: str
) -> None:
    """Checks that a nitro group written with its usual [O-] is still
    neutralized and handed to the Nitro site, unlike an N-oxide."""

    check_protonation("C[N+](=O)[O-]", ph, [expected], [state])


@pytest.mark.parametrize(
    "smiles", ["[O-][N+](=O)[O-]", "O[N+](=O)[O-]"], ids=["nitrate", "nitric_acid"]
)
@pytest.mark.parametrize(
    "ph, expected, state",
    [
        (VERY_ACIDIC_PH, "O[N+](=O)[O-]", "PROTONATED"),
        (7.0, "[O-][N+](=O)[O-]", "DEPROTONATED"),
        (VERY_BASIC_PH, "[O-][N+](=O)[O-]", "DEPROTONATED"),
    ],
    ids=["very_acidic", "neutral", "very_basic"],
)
def test_nitrate_loses_both_protons(
    smiles: str, ph: float, expected: str, state: str
) -> None:
    """Checks that nitrate is not reported as nitric acid. Neutralization
    protonated both of its O-, and the Nitro pattern then restored only one,
    because the other OH was locked as context. Nitrate counterions came out
    as neutral HNO3 at pH 7."""

    check_protonation(smiles, ph, [expected], [state])


@pytest.mark.parametrize(
    "ph, expected, state",
    [
        (VERY_ACIDIC_PH, "CO[N+](=O)O", "PROTONATED"),
        (VERY_BASIC_PH, "CO[N+](=O)[O-]", "DEPROTONATED"),
    ],
    ids=["very_acidic", "very_basic"],
)
def test_nitrate_ester_still_handled(ph: float, expected: str, state: str) -> None:
    """Checks that the nitrate exception in neutralization does not reach
    nitrate esters, whose single O- must still be protonated for Nitro."""

    check_protonation("CO[N+](=O)[O-]", ph, [expected], [state])


@pytest.mark.parametrize(
    "charged, neutral",
    [
        ("C[S+](C)[O-]", "CS(C)=O"),
        ("C[S+2](C)([O-])[O-]", "CS(C)(=O)=O"),
        ("C[P+](C)(C)[O-]", "CP(C)(C)=O"),
        ("CO[P+]([O-])([O-])[O-]", "COP(=O)(O)O"),
    ],
    ids=["sulfoxide", "sulfone", "phosphine_oxide", "phosphate"],
)
def test_charge_separated_oxides_are_collapsed(
    charged: str, neutral: str, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that S+-O- and P+-O- are restored to S=O and P=O. The O- rule
    protonated these oxygens instead, so DMSO came out as the hydroxysulfonium
    cation C[S+](C)O at every pH."""

    record = dimorphite_dl.LoadSMIFile(StringIO(charged)).next()
    assert record["smiles"] == canonical_smiles(neutral)


CHARGE_ON_NON_H_NITROGEN = [
    # [charged input, neutral form, id]
    ["C[n+]1cc[nH]c1", "Cn1ccnc1", "imidazolium"],
    ["C[n+]1ccc[nH]1", "Cn1cccn1", "pyrazolium"],
    ["CC(N)=[N+](C)C", "CC(=N)N(C)C", "amidinium"],
    ["CN(C)C(N)=[N+](C)C", "CN(C)C(=N)N(C)C", "guanidinium"],
]


@pytest.mark.parametrize(
    "charged, neutral",
    [group[:2] for group in CHARGE_ON_NON_H_NITROGEN],
    ids=[group[2] for group in CHARGE_ON_NON_H_NITROGEN],
)
def test_charge_on_non_h_nitrogen_is_neutralized(
    charged: str, neutral: str, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that a cation drawn with the charge on a nitrogen that has no H
    is neutralized. Only N+-H lost its proton, so these resonance forms kept
    their charge and were never matched by the neutral site patterns."""

    record = dimorphite_dl.LoadSMIFile(StringIO(charged)).next()
    assert record["smiles"] == canonical_smiles(neutral)


@pytest.mark.parametrize(
    "charged, neutral",
    [group[:2] for group in CHARGE_ON_NON_H_NITROGEN],
    ids=[group[2] for group in CHARGE_ON_NON_H_NITROGEN],
)
@pytest.mark.parametrize(
    "ph", [VERY_ACIDIC_PH, VERY_BASIC_PH], ids=["very_acidic", "very_basic"]
)
def test_charge_on_non_h_nitrogen_matches_neutral_output(
    charged: str, neutral: str, ph: float
) -> None:
    """Checks that the resonance form of the input does not change the
    output. The charged imidazolium form was emitted at every pH, and at very
    basic pH its N-H was even removed to give C[n+]1cc[n-]c1."""

    expected = sorted(
        [normalize_smiles(line[0])] + line[1:]
        for line in protonate(neutral, ph, DEFAULT_PKA_PRECISION)
    )
    output = sorted(
        [normalize_smiles(line[0])] + line[1:]
        for line in protonate(charged, ph, DEFAULT_PKA_PRECISION)
    )
    assert output == expected


@pytest.mark.parametrize(
    "smiles",
    ["C[n+]1ccn(C)c1", "CC(N(C)C)=[N+](C)C"],
    ids=["dimethylimidazolium", "peralkyl_amidinium"],
)
def test_fully_substituted_cation_stays_charged(
    smiles: str, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that the charge-shift rules leave a cation with no N-H alone,
    since it has no proton to lose and really is permanently charged."""

    record = dimorphite_dl.LoadSMIFile(StringIO(smiles)).next()
    assert record["smiles"] == canonical_smiles(smiles)


DIAZO_COMPOUNDS = [
    # [smiles, id]
    ["CCOC(=O)C=[N+]=[N-]", "ethyl_diazoacetate"],
    ["[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)O", "azaserine"],
    ["[N-]=[N+]=CC(=O)CC[C@H](N)C(=O)O", "DON"],
]


@pytest.mark.parametrize(
    "smiles",
    [group[0] for group in DIAZO_COMPOUNDS],
    ids=[group[1] for group in DIAZO_COMPOUNDS],
)
def test_diazo_input_is_not_neutralized(
    smiles: str, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that the terminal N- of a diazo group is left alone. The rule
    that turns azide N- into N-H also matched it, giving C=[N+]=N, which no
    site pattern deprotonates."""

    record = dimorphite_dl.LoadSMIFile(StringIO(smiles)).next()
    assert record["smiles"] == canonical_smiles(smiles)


@pytest.mark.parametrize(
    "ph", [VERY_ACIDIC_PH, 7.4, VERY_BASIC_PH], ids=["very_acidic", "7.4", "very_basic"]
)
def test_diazo_without_sites_is_unchanged(
    ph: float, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that a diazo compound with no ionizable site comes out neutral
    at every pH, rather than as a permanent +1 cation."""

    smiles = "CCOC(=O)C=[N+]=[N-]"
    output = protonate(smiles, ph, DEFAULT_PKA_PRECISION)
    assert [line[0] for line in output] == [canonical_smiles(smiles)], output


@pytest.mark.parametrize(
    "ph, expected, state",
    [
        (VERY_ACIDIC_PH, "[N-]=[N+]=CC(=O)OC[C@H]([NH3+])C(=O)O", "PROTONATED"),
        (VERY_BASIC_PH, "[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)[O-]", "DEPROTONATED"),
    ],
    ids=["very_acidic", "very_basic"],
)
def test_diazo_does_not_shift_net_charge_of_other_sites(
    ph: float, expected: str, state: str
) -> None:
    """Checks that azaserine's amine and acid are protonated as usual while
    its diazo group stays neutral, so the net charge is +1 and -1 rather than
    +2 and 0."""

    check_protonation("[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)O", ph, [expected], [state])


THIOLATE_GROUPS = [
    # [input smiles, protonated, deprotonated, id]
    ["CC(C)(C)[S-]", "CC(C)(C)S", "CC(C)(C)[S-]", "Thiol"],
    ["[S-]c1ccccc1", "Sc1ccccc1", "[S-]c1ccccc1", "Phenyl_Thiol"],
    ["CC(=O)[S-]", "CC(=O)S", "CC(=O)[S-]", "Thioic_acid"],
]


@pytest.mark.parametrize("group", THIOLATE_GROUPS, ids=group_id)
@pytest.mark.parametrize(
    "ph, state",
    [(VERY_ACIDIC_PH, "PROTONATED"), (VERY_BASIC_PH, "DEPROTONATED")],
    ids=["very_acidic", "very_basic"],
)
def test_thiolate_input_is_neutralized(group: List[str], ph: float, state: str) -> None:
    """Checks that sulfur given already deprotonated is still recognized as a
    site. The thiol patterns require S-H, so [S-] inputs kept their charge at
    any pH."""

    smiles, protonated, deprotonated, _ = group
    expected = protonated if state == "PROTONATED" else deprotonated
    check_protonation(smiles, ph, [expected], [state])


def test_multiple_both_sites_enumerate_every_combination(
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that two BOTH sites combine into all four states rather than
    only the fully protonated and fully deprotonated ones."""

    output = protonate("NCCCC(=O)O", 7.0, EVERY_SITE_BOTH_PRECISION)

    expected = [
        "NCCCC(=O)O",
        "[NH3+]CCCC(=O)O",
        "NCCCC(=O)[O-]",
        "[NH3+]CCCC(=O)[O-]",
    ]
    assert sorted(line[0] for line in output) == sorted(
        canonical_smiles(s) for s in expected
    ), output
    assert all(line[1:] == ["BOTH", "BOTH"] for line in output), output


def test_unbuildable_state_is_reported(capsys: pytest.CaptureFixture[str]) -> None:
    """Checks that a state with no valid structure is reported. The parent
    is emitted in its place under the requested label, so the pyridone below
    came out as one BOTH line with nothing saying its protonated form had
    failed."""

    output = protonate("O=c1cc[nH]cc1", 7.0, EVERY_SITE_BOTH_PRECISION)

    assert [line[1:] for line in output] == [["BOTH"]], output
    err = capsys.readouterr().err
    assert err.count("No valid protonated state") == 1, err
    assert "Aromatic_nitrogen_protonated" in err, err


def test_buildable_states_are_not_reported(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Checks that the warning is specific to failed states, so it cannot
    become noise on ordinary input."""

    protonate("NCCCC(=O)O", 7.0, EVERY_SITE_BOTH_PRECISION)

    assert "No valid" not in capsys.readouterr().err


def test_many_both_sites_uncapped(capsys: pytest.CaptureFixture[str]) -> None:
    """Checks that every combination is produced when the cap is not hit, and
    that no truncation warning is printed."""

    output = protonate(
        EIGHT_CARBOXYLS,
        7.0,
        EVERY_SITE_BOTH_PRECISION,
        max_variants=EIGHT_CARBOXYL_STATES,
    )

    output_smiles = [line[0] for line in output]
    assert len(output_smiles) == EIGHT_CARBOXYL_STATES
    assert len(set(output_smiles)) == EIGHT_CARBOXYL_STATES
    assert all(Chem.MolFromSmiles(s) is not None for s in output_smiles)
    assert "Limited number of variants" not in capsys.readouterr().err


@pytest.mark.parametrize(
    "max_variants",
    [1, 10, dimorphite_dl.DEFAULT_MAX_VARIANTS],
    ids=["one", "ten", "default"],
)
def test_many_both_sites_capped(
    max_variants: int,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Checks that the cap limits the output, that the working list never
    grows past twice the cap (so memory stays bounded), and that the
    truncation is reported."""

    original = dimorphite_dl.ProtSubstructFuncs.protonate_site
    input_sizes: List[int] = []

    def recording_protonate_site(
        mols: List[Chem.Mol], site: Tuple[int, str, str, float]
    ) -> List[Chem.Mol]:
        """Records how many states each site is applied to.

        Args:
            mols: The states going into this site.
            site: The site being protonated.

        Returns:
            The original function's result.
        """

        input_sizes.append(len(mols))
        return original(mols, site)

    monkeypatch.setattr(
        dimorphite_dl.ProtSubstructFuncs,
        "protonate_site",
        staticmethod(recording_protonate_site),
    )

    output = protonate(
        EIGHT_CARBOXYLS, 7.0, EVERY_SITE_BOTH_PRECISION, max_variants=max_variants
    )

    output_smiles = [line[0] for line in output]
    assert len(output_smiles) == max_variants
    assert len(set(output_smiles)) == max_variants
    assert all(Chem.MolFromSmiles(s) is not None for s in output_smiles)

    # Every site is still assigned after the cap is reached.
    assert len(input_sizes) == 8
    assert max(input_sizes) <= max_variants

    assert "Limited number of variants" in capsys.readouterr().err


def load_input_mol(smiles: str) -> Chem.Mol:
    """Reads one SMILES string the way Protonate does (neutralized and
    canonicalized), so site tests see the molecule the sites are found on.

    Args:
        smiles: The input SMILES string.

    Returns:
        The Mol that Protonate would search for sites.
    """

    smi = dimorphite_dl.LoadSMIFile(StringIO(smiles)).next()["smiles"]
    return dimorphite_dl.UtilFuncs.convert_smiles_str_to_mol(smi)


@pytest.mark.parametrize(
    "smiles",
    [group[0] for group in SINGLE_SITE_GROUPS + PHOSPHORUS_GROUPS],
    ids=[group_id(group) for group in SINGLE_SITE_GROUPS + PHOSPHORUS_GROUPS],
)
def test_every_site_is_an_atom_of_the_input_mol(smiles: str) -> None:
    """Checks that no shipped pattern puts a site on a hydrogen that AddHs
    created. Sites are found on the H-added copy but charged on the Mol
    without those Hs, so such a site would name an atom that is not there."""

    mol = load_input_mol(smiles)
    subs = (
        dimorphite_dl.ProtSubstructFuncs.load_protonation_substructs_calc_state_for_ph()
    )
    sites = dimorphite_dl.ProtSubstructFuncs.get_prot_sites_and_target_states(mol, subs)

    assert sites != [], smiles
    for idx, _, name, _ in sites:
        assert idx < mol.GetNumAtoms(), (name, idx)
        assert mol.GetAtomWithIdx(idx).GetAtomicNum() != 1, (name, idx)


def test_site_on_added_hydrogen_raises() -> None:
    """Checks that a pattern whose site is an AddHs hydrogen raises instead
    of charging whatever atom happens to have that index."""

    smarts = "[OX2]-[#1]"
    subs: List[dimorphite_dl.SiteSubstruct] = [
        {
            "name": "Hydroxyl_hydrogen",
            "smart": smarts,
            "mol": Chem.MolFromSmarts(smarts),
            "prot_states_for_pH": [("1", "BOTH", 7.0)],
        }
    ]

    with pytest.raises(ValueError, match="Hydroxyl_hydrogen"):
        dimorphite_dl.ProtSubstructFuncs.get_prot_sites_and_target_states(
            Chem.MolFromSmiles("CCO"), subs
        )


@pytest.mark.parametrize(
    "group",
    SINGLE_SITE_GROUPS + PHOSPHORUS_GROUPS,
    ids=[group_id(group) for group in SINGLE_SITE_GROUPS + PHOSPHORUS_GROUPS],
)
def test_output_does_not_depend_on_input_charge_state(group: List[str]) -> None:
    """Checks that each group gives the same states and labels whether it is
    written neutral, protonated, or deprotonated. neutralize_mol handles only
    the charged forms it has a rule for, so a missing rule shows up here as a
    charged input that keeps its charge or loses its site."""

    smiles = group[0]
    expected = sorted(
        [normalize_smiles(line[0])] + line[1:]
        for line in protonate(smiles, 7.0, EVERY_SITE_BOTH_PRECISION)
    )

    for charged_form in group[1:-1]:
        output = sorted(
            [normalize_smiles(line[0])] + line[1:]
            for line in protonate(charged_form, 7.0, EVERY_SITE_BOTH_PRECISION)
        )
        assert output == expected, charged_form


def formal_charge(smiles: str) -> int:
    """Sums a molecule's formal charges, which for the carboxyl tests is minus
    the number of deprotonated groups.

    Args:
        smiles: A valid SMILES string.

    Returns:
        The molecule's net formal charge.
    """

    return sum(atom.GetFormalCharge() for atom in Chem.MolFromSmiles(smiles).GetAtoms())


@pytest.mark.parametrize(
    "max_variants, expected",
    [
        (1, ["[NH3+]CCCC(=O)[O-]"]),
        (2, ["[NH3+]CCCC(=O)[O-]", "NCCCC(=O)[O-]"]),
        (3, ["[NH3+]CCCC(=O)[O-]", "NCCCC(=O)[O-]", "[NH3+]CCCC(=O)O"]),
    ],
    ids=["one", "two", "three"],
)
def test_max_variants_keeps_most_probable_states(
    max_variants: int, expected: List[str], canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that truncation keeps the likeliest states at the pH. At pH 7
    the amine (pKa about 8.2) is mostly protonated and the carboxyl (about
    3.5) mostly deprotonated. Taking the first states generated instead kept
    the neutral amine, because deprotonated copies are generated first."""

    output = protonate("NCCCC(=O)O", 7.0, EVERY_SITE_BOTH_PRECISION, max_variants)

    assert sorted(line[0] for line in output) == sorted(
        canonical_smiles(s) for s in expected
    ), output


def test_max_variants_does_not_pin_later_sites() -> None:
    """Checks that with eight equivalent carboxyls at pH 7 and room for nine
    states, the cap keeps the fully deprotonated state and all eight singly
    protonated ones. Taking the first states generated instead pinned every
    carboxyl after the cap was first reached to its deprotonated form."""

    output = protonate(EIGHT_CARBOXYLS, 7.0, EVERY_SITE_BOTH_PRECISION, 9)

    charges = sorted(formal_charge(line[0]) for line in output)
    assert charges == [-8] + [-7] * 8, output


@pytest.mark.parametrize(
    "pka, ph", [(7.0, 7.0), (4.0, 7.0), (10.0, 7.0), (-1000.0, 7.4)]
)
def test_charge_log_probabilities_are_fractions(pka: float, ph: float) -> None:
    """Checks that a BOTH site's deprotonated and protonated fractions sum to
    one, split evenly at the pKa, and do not overflow at Nitro's pKa."""

    deprotonated, protonated = (
        dimorphite_dl.ProtSubstructFuncs.charge_log_probabilities("BOTH", pka, ph)
    )

    assert 10**deprotonated + 10**protonated == pytest.approx(1.0)
    if pka == ph:
        assert deprotonated == pytest.approx(protonated)
    else:
        assert (protonated > deprotonated) == (pka > ph)


@pytest.mark.parametrize("state", ["PROTONATED", "DEPROTONATED"])
def test_single_state_site_does_not_affect_ranking(state: str) -> None:
    """Checks that a site with only one charge scores zero, since it shifts
    every state equally."""

    assert dimorphite_dl.ProtSubstructFuncs.charge_log_probabilities(
        state, 8.0, 7.0
    ) == [0.0]


@pytest.mark.parametrize(
    "mean, std, min_ph, max_ph",
    [(7.4, 0.1, 8.4, 6.4), (7.4, -0.1, 6.4, 8.4)],
    ids=["inverted_ph_range", "negative_std"],
)
def test_define_protonation_state_rejects_inverted_intervals(
    mean: float, std: float, min_ph: float, max_ph: float
) -> None:
    """Checks that an inverted pH or pKa interval raises instead of silently
    returning the wrong state."""

    with pytest.raises(ValueError):
        dimorphite_dl.ProtSubstructFuncs.define_protonation_state(
            mean, std, min_ph, max_ph
        )


@pytest.mark.parametrize(
    "params",
    [
        {"min_ph": 8.4, "max_ph": 6.4},
        {"pka_precision": -1.0},
        {"max_variants": 0},
        {"max_ph": float("nan")},
        {"min_ph": float("-inf")},
        {"pka_precision": float("nan")},
    ],
    ids=[
        "inverted_ph_range",
        "negative_precision",
        "zero_max_variants",
        "nan_max_ph",
        "infinite_min_ph",
        "nan_precision",
    ],
)
def test_protonate_rejects_invalid_ranges(params: Dict[str, float]) -> None:
    """Checks that bad user parameters are rejected before any protonation."""

    args: Dict[str, Union[str, float]] = {"smiles": "CCCN"}
    args.update(params)
    with pytest.raises(ValueError):
        dimorphite_dl.Protonate(args)
