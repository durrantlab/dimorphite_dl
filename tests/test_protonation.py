"""Checks the protonation states Dimorphite-DL assigns to each ionizable
group, alone and in combination."""

import os
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

    expected = set(normalize_smiles(s) for s in expected_smiles)
    assert len(output) == len(expected), output
    assert set(normalize_smiles(s) for s in output_smiles) <= expected, output
    assert set(label for line in output for label in line[1:]) <= set(labels), output
    return output


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
        mols: List[Chem.Mol], site: Tuple[int, str, str]
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
    [{"min_ph": 8.4, "max_ph": 6.4}, {"pka_precision": -1.0}, {"max_variants": 0}],
    ids=["inverted_ph_range", "negative_precision", "zero_max_variants"],
)
def test_protonate_rejects_invalid_ranges(params: Dict[str, float]) -> None:
    """Checks that bad user parameters are rejected before any protonation."""

    args: Dict[str, Union[str, float]] = {"smiles": "CCCN"}
    args.update(params)
    with pytest.raises(ValueError):
        dimorphite_dl.Protonate(args)
