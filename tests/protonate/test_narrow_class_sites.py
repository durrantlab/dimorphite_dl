"""Checks the narrow site entries added alongside the broad trained ones:
the alkyl, benzylic, beta-hydroxy and fluoroalkyl amines, N-alkyl anilines,
unactivated phenols, alkanethiols, primary sulfonamides, C-substituted
amidines, and N-nitro amides.

Each of those entries exists to give one state where its parent rule gave two,
so one half of this file checks that the collapse happens. The other half is
the more important one: every entry carries exclusions that hand specific
compounds back to the parent rule, and those compounds are correct today only
because of the exclusions. Nothing else records that, so editing a pattern
could silently reclaim them.
"""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles

DEFAULT_PH_MIN = 6.4
DEFAULT_PH_MAX = 8.4


def canonical_set(smiles_list: list[str]) -> set[str]:
    """Compare outputs by structure rather than by SMILES spelling.

    Args:
        smiles_list: SMILES strings to canonicalize.

    Returns:
        The set of canonical SMILES.
    """
    return {Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in smiles_list}


def protonate_default(smiles: str) -> list[str]:
    """Runs the default pH window, which is what the narrow entries were sized
    against.

    Args:
        smiles: The input SMILES string.

    Returns:
        The output SMILES.
    """
    return protonate_smiles(smiles, ph_min=DEFAULT_PH_MIN, ph_max=DEFAULT_PH_MAX)


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        pytest.param("CN", "C[NH3+]", id="methylamine_10.62"),
        pytest.param("C1CCNCC1", "C1CC[NH2+]CC1", id="piperidine_11.22"),
        pytest.param("CCN(CC)CC", "CC[NH+](CC)CC", id="triethylamine_10.65"),
        pytest.param("C1CN2CCC1CC2", "C1C[NH+]2CCC1CC2", id="quinuclidine_11.0"),
        pytest.param("NCc1ccccc1", "[NH3+]Cc1ccccc1", id="benzylamine_9.34"),
        pytest.param("NCCO", "[NH3+]CCO", id="ethanolamine_9.50"),
        pytest.param(
            "CC(C)NCC(O)COc1cccc2ccccc12",
            "CC(C)[NH2+]CC(O)COc1cccc2ccccc12",
            id="propranolol_9.42",
        ),
        pytest.param("CC(N)=N", "CC(N)=[NH2+]", id="acetamidine_12.52"),
        pytest.param("NC(=N)c1ccccc1", "NC(=[NH2+])c1ccccc1", id="benzamidine_11.6"),
    ],
)
def test_base_is_cation_only(smiles: str, expected: str) -> None:
    """A base well above ph_max should give the cation alone. Each of these
    was returned in both states before its entry existed, because the parent
    rule's standard deviation spans the whole default window.

    Args:
        smiles: The neutral input.
        expected: The only state that should come back.
    """
    output = protonate_default(smiles)
    assert canonical_set(output) == canonical_set([expected]), output


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        pytest.param("c1ccccc1O", "Oc1ccccc1", id="phenol_9.99"),
        pytest.param("CCS", "CCS", id="ethanethiol_10.50"),
        pytest.param(
            "NS(=O)(=O)c1ccccc1", "NS(=O)(=O)c1ccccc1", id="benzenesulfonamide_10.1"
        ),
        pytest.param("CS(=O)(=O)N", "CS(N)(=O)=O", id="methanesulfonamide_10.8"),
        pytest.param("CNc1ccccc1", "CNc1ccccc1", id="N_methylaniline_4.84"),
        pytest.param("FC(F)(F)CN", "NCC(F)(F)F", id="trifluoroethylamine_5.7"),
    ],
)
def test_weak_site_is_neutral_only(smiles: str, expected: str) -> None:
    """A weak acid above the window, or a weak base below it, should give the
    neutral molecule alone.

    Args:
        smiles: The input.
        expected: The only state that should come back.
    """
    output = protonate_default(smiles)
    assert canonical_set(output) == canonical_set([expected]), output


def test_n_nitro_amide_is_anion_and_keeps_its_nitro_group() -> None:
    """N-nitrourethane has a pKa of 3.28, so the N-H should be gone. It came
    back neutral until Nitro stopped claiming the atom its nitro group hangs
    from, which had protected this nitrogen. Both sites have to fire on the
    one molecule, so the nitro group must also survive as the anion."""
    output = protonate_default("CCOC(=O)N([H])[N+](=O)[O-]")
    assert canonical_set(output) == canonical_set(["CCOC(=O)[N-][N+](=O)[O-]"]), output


def test_nitroguanidine_is_not_basic() -> None:
    """A nitro group on a guanidine drops the conjugate acid from 13.6 to
    -0.98, so the neutral molecule is the only state at physiological pH.
    AmidineGuanidine1 matches any carbon with three nitrogens and was
    returning the cation. Nitroarginine is the case that surfaced it."""
    output = protonate_default("CCCNC(=N)N[N+](=O)[O-]")
    assert canonical_set(output) == canonical_set(["CCCNC(=N)N[N+](=O)[O-]"]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param(r"O=[N+]([O-])/N=C1\NCOCN1Cc1cnc(Cl)s1", id="thiamethoxam"),
        pytest.param("O=[N+]([O-])N=C1NCCN1Cc1ccc(Cl)nc1", id="imidacloprid"),
    ],
)
def test_nitroimine_guanidine_is_not_basic(smiles: str) -> None:
    """The neonicotinoids carry the nitro group on the imine nitrogen, where
    it costs about twelve units of basicity (imidacloprid's conjugate acid is
    1.56). They fell through to AmidineGuanidine2 and were returned in both
    states, which only became visible once the nitro entry stopped protecting
    the atoms around it.

    The total number of states is not checked, because both compounds also
    carry a chloro-heteroaryl ring whose nitrogen is a site in its own right;
    only the nitroimine nitrogen is.

    Args:
        smiles: The neutral input.
    """
    neutral = Chem.MolFromSmarts(
        "[NX2+0;H0;$(N-[N+](=[OX1])[OX1-]),$(N-[NX3](=[OX1])=[OX1])]=[#6]"
    )
    assert neutral is not None
    for out in protonate_default(smiles):
        mol = Chem.MolFromSmiles(out)
        assert mol is not None, out
        assert mol.HasSubstructMatch(neutral), out


def test_plain_guanidine_is_still_a_cation() -> None:
    """Nitroguanidine sits ahead of AmidineGuanidine1, so it must not poach
    guanidines that carry no nitro group."""
    output = protonate_default("NC(N)=N")
    assert canonical_set(output) == canonical_set(["NC(N)=[NH2+]"]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CCCNC(=N)N[N+](=O)[O-]", id="nitroarginine_fragment"),
        pytest.param("C[N+](=O)[O-]", id="nitromethane"),
        pytest.param("CC(C)[N+](=O)[O-]", id="2-nitropropane"),
        pytest.param("Oc1ccc(cc1)[N+](=O)[O-]", id="4-nitrophenol"),
    ],
)
def test_nitro_group_survives_as_the_anion(smiles: str) -> None:
    """Nitro now names its substituent in a recursive constraint rather than
    matching it. The group itself must still be deprotonated, since the
    neutralizer protonates it on the way in.

    Args:
        smiles: A molecule with a charge-separated nitro group.
    """
    # The alternatives have to be recursive primitives inside one atom
    # expression; a comma between two bonded expressions is not valid SMARTS
    # and MolFromSmarts returns None for it.
    nitro = Chem.MolFromSmarts("[$([N+](=[OX1])[OX1-]),$([NX3](=[OX1])=[OX1])]")
    assert nitro is not None
    for out in protonate_default(smiles):
        mol = Chem.MolFromSmiles(out)
        assert mol is not None, out
        assert mol.HasSubstructMatch(nitro), out


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("OC(=O)C(N)CS", id="cysteine_8.3"),
        pytest.param("C1COCCN1", id="morpholine_8.36"),
        pytest.param("C1CN1", id="aziridine_7.98"),
        pytest.param("NC(CO)(CO)CO", id="tris_8.10"),
        pytest.param("OCCNCCO", id="diethanolamine_9.00"),
        pytest.param("OCCN(CCO)CCO", id="triethanolamine_7.77"),
        pytest.param("CCN(CC)CC(=O)Nc1c(C)cccc1C", id="lidocaine_7.86"),
        pytest.param("CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1", id="fentanyl_8.44"),
        pytest.param(
            "CCC(=O)N(c1ccccc1)C1(COC)CCN(CCc2cccs2)CC1", id="sufentanil_8.01"
        ),
        pytest.param(
            "CCOC(=O)C1CN(Cc2ccccc2)CC1C(F)(F)F",
            id="cf3_benzylpyrrolidine_ester_5.40",
        ),
        pytest.param("NCCN", id="ethylenediamine_9.98"),
        pytest.param("C1CN2CCN1CC2", id="dabco_8.19"),
        pytest.param("C=CCN(CC=C)CC=C", id="triallylamine_8.31"),
        pytest.param("C#CCN", id="propargylamine_7.05"),
        pytest.param("Oc1ccc(cc1)C(F)(F)F", id="4-trifluoromethylphenol_8.68"),
        pytest.param("CC(=O)c1ccc(O)cc1", id="4-hydroxyacetophenone_8.0"),
        pytest.param("Oc1ccc2ccccc2c1", id="2-naphthol_9.51"),
        pytest.param(
            "NC(=N)NC(C)=O",
            id="acetylguanidine_8.33",
            marks=pytest.mark.xfail(
                strict=False,
                reason=(
                    "Pre-existing: acylation drops guanidine from 13.71 to 8.33, "
                    "but the carbon carries three nitrogens so AmidineGuanidine1 "
                    "claims it at 12.03 and returns the cation alone. The case is "
                    "kept because it also shows Amidine_carbon_substituted, which "
                    "requires carbon or hydrogen as the third substituent, is not "
                    "poaching it. Fixing it means carving acylguanidines out of a "
                    "trained entry."
                ),
            ),
        ),
    ],
)
def test_excluded_compound_keeps_both_states(smiles: str) -> None:
    """Each of these sits inside the default window, so both states are right,
    and each is excluded from a narrow entry on a specific ground: a beta or
    gamma heteroatom, a three-membered ring, a second hydroxyl, an unsaturated
    or activated neighbor, a fused or substituted ring, or an acyl group. If a
    narrow entry reclaims one, this test is how you find out.

    Args:
        smiles: A compound deliberately left on a trained rule.
    """
    output = protonate_default(smiles)
    assert len(output) > 1, output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("Oc1ccc(cc1)[N+](=O)[O-]", id="4-nitrophenol_7.14"),
        pytest.param("N#Cc1ccc(O)cc1", id="4-cyanophenol_7.95"),
        pytest.param("Oc1ccccc1O", id="catechol_9.25"),
        pytest.param("NS(=O)(=O)c1ccc(N)cc1", id="sulfanilamide"),
    ],
)
def test_ring_test_keeps_activated_aromatics_on_the_generic_rule(
    smiles: str,
) -> None:
    """Phenol_unactivated and Primary_sulfonamide both require a ring whose
    other positions carry only hydrogen or a carbon bearing nothing but carbon
    and hydrogen. These rings all fail that test and must reach the broad
    entry, whose wide window is what makes them come out in both states.

    Args:
        smiles: An aromatic bearing a heteroatom or activated substituent.
    """
    output = protonate_default(smiles)
    assert len(output) > 1, output


@pytest.mark.parametrize(
    ("neutral", "charged"),
    [
        pytest.param("CN", "C[NH3+]", id="methylamine"),
        pytest.param("C1CCNCC1", "C1CC[NH2+]CC1", id="piperidine"),
        pytest.param("NCCO", "[NH3+]CCO", id="ethanolamine"),
        pytest.param("NCc1ccccc1", "[NH3+]Cc1ccccc1", id="benzylamine"),
        pytest.param("CC(N)=N", "CC(N)=[NH2+]", id="acetamidine"),
        pytest.param("c1ccccc1O", "[O-]c1ccccc1", id="phenol"),
        pytest.param("CCS", "CC[S-]", id="ethanethiol"),
        pytest.param("NS(=O)(=O)c1ccccc1", "[NH-]S(=O)(=O)c1ccccc1", id="sulfonamide"),
    ],
)
def test_answer_does_not_depend_on_input_charge(neutral: str, charged: str) -> None:
    """Several of the narrow patterns specify a formal charge of zero or an
    explicit hydrogen count, so a pre-charged input could miss them and fall
    through to the parent rule. Both spellings must give the same answer.

    Args:
        neutral: The neutral form of the molecule.
        charged: The same molecule drawn in its ionized form.
    """
    assert canonical_set(protonate_default(neutral)) == canonical_set(
        protonate_default(charged)
    )
