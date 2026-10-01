import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from dimorphite_dl import protonate_smiles
from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.neutralize import MoleculeNeutralizer

VERY_ACIDIC_PH = -10000000.0
VERY_BASIC_PH = 10000000.0


@pytest.mark.parametrize(
    ("input_smiles", "exp_azides", "exp_neutral", "exp_canonical"),
    [
        ("C#CCO", "C#CCO", "C#CCO", "C#CCO"),
        ("Brc1cc[nH+]cc1", "Brc1cc[nH+]cc1", "Brc1ccncc1", "Brc1ccncc1"),
        ("C-N=[N+]=[N@H]", "C-N=[N+]=[N@H]", "CN=[N+]=N", "CN=[N+]=N"),
        ("O=P(O)(O)OCCCC", "O=P(O)(O)OCCCC", "CCCCOP(=O)(O)O", "CCCCOP(=O)(O)O"),
    ],
)
def test_molecule_preparation_steps(
    input_smiles, exp_azides, exp_neutral, exp_canonical
):
    mol = MoleculeRecord(input_smiles)
    mol.process_azides()
    assert (
        mol.smiles == exp_azides
    ), f"after process_azides: got {mol.smiles!r}, expected {exp_azides!r}"
    mol.neutralize()
    assert (
        mol.smiles == exp_neutral
    ), f"after neutralize: got {mol.smiles!r}, expected {exp_neutral!r}"
    mol.make_canonical()
    assert (
        mol.smiles == exp_canonical
    ), f"after make_canonical: got {mol.smiles!r}, expected {exp_canonical!r}"


def canonical(smiles: str) -> str:
    """Canonicalizes so expected values can be written in any valid form.

    Args:
        smiles: A valid SMILES string.

    Returns:
        RDKit's canonical isomeric SMILES for smiles.
    """

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


def neutralized(smiles: str) -> str:
    """Runs only the neutralization step, the way protonate_smiles prepares
    each input.

    Args:
        smiles: The input SMILES string.

    Returns:
        The canonical SMILES after neutralization.
    """

    record = MoleculeRecord(smiles)
    record.process_azides()
    record.neutralize()
    return canonical(record.smiles)


def protonated(smiles: str, ph: float) -> list[str]:
    """Runs the full pipeline at a single pH.

    Args:
        smiles: The input SMILES string.
        ph: Used as both the minimum and maximum pH.

    Returns:
        The sorted canonical SMILES of every output state.
    """

    output = protonate_smiles(smiles, ph_min=ph, ph_max=ph, precision=0.5)
    return sorted(canonical(line) for line in output)


@pytest.mark.parametrize(
    ("charged", "neutral"),
    [
        ("C[S+](C)[O-]", "CS(C)=O"),
        ("C[S+2](C)([O-])[O-]", "CS(C)(=O)=O"),
        ("C[P+](C)(C)[O-]", "CP(C)(C)=O"),
        ("CO[P+]([O-])([O-])[O-]", "COP(=O)(O)O"),
    ],
    ids=["sulfoxide", "sulfone", "phosphine_oxide", "phosphate"],
)
def test_charge_separated_oxides_are_collapsed(charged: str, neutral: str) -> None:
    """Checks that S+-O- and P+-O- are restored to S=O and P=O. The O- rule
    protonated these oxygens instead, so DMSO came out as the hydroxysulfonium
    cation C[S+](C)O at every pH."""

    assert neutralized(charged) == canonical(neutral)


@pytest.mark.parametrize(
    "smiles",
    ["C[N+](C)(C)[O-]", "[O-][n+]1ccccc1", "CC=[N+](C)[O-]"],
    ids=["amine_oxide", "pyridine_n_oxide", "nitrone"],
)
def test_n_oxide_keeps_oxide(smiles: str) -> None:
    """Checks that an N-oxide is not turned into an N-hydroxy cation. No site
    restores the O-, so amine oxides came out as [N+]-OH at every pH."""

    assert neutralized(smiles) == canonical(smiles)
    assert protonated(smiles, 7.0) == [canonical(smiles)]


@pytest.mark.parametrize(
    ("ph", "expected"),
    [(VERY_ACIDIC_PH, "C[N+](=O)O"), (VERY_BASIC_PH, "C[N+](=O)[O-]")],
    ids=["very_acidic", "very_basic"],
)
def test_charged_nitro_input_still_handled(ph: float, expected: str) -> None:
    """Checks that the N-oxide exception does not reach nitro groups, whose
    O- must still be protonated so the Nitro site can match."""

    assert protonated("C[N+](=O)[O-]", ph) == [canonical(expected)]


@pytest.mark.parametrize(
    "smiles", ["[O-][N+](=O)[O-]", "O[N+](=O)[O-]"], ids=["nitrate", "nitric_acid"]
)
@pytest.mark.parametrize(
    ("ph", "expected"),
    [
        (VERY_ACIDIC_PH, "O[N+](=O)[O-]"),
        (7.0, "[O-][N+](=O)[O-]"),
        (VERY_BASIC_PH, "[O-][N+](=O)[O-]"),
    ],
    ids=["very_acidic", "neutral", "very_basic"],
)
def test_nitrate_is_not_reported_as_nitric_acid(
    smiles: str, ph: float, expected: str
) -> None:
    """Checks that nitrate keeps a single ionizable OH. Neutralization
    protonated both O-, giving O[N+](=O)O, which is not nitric acid and has
    two Nitro sites instead of one."""

    assert protonated(smiles, ph) == [canonical(expected)]


@pytest.mark.parametrize(
    ("ph", "expected"),
    [(VERY_ACIDIC_PH, "CO[N+](=O)O"), (VERY_BASIC_PH, "CO[N+](=O)[O-]")],
    ids=["very_acidic", "very_basic"],
)
def test_nitrate_ester_still_handled(ph: float, expected: str) -> None:
    """Checks that the nitrate exception does not reach nitrate esters, whose
    single O- must still be protonated for Nitro."""

    assert protonated("CO[N+](=O)[O-]", ph) == [canonical(expected)]


@pytest.mark.parametrize(
    ("smiles", "protonated_form", "deprotonated_form"),
    [
        ("CC(C)(C)[S-]", "CC(C)(C)S", "CC(C)(C)[S-]"),
        ("[S-]c1ccccc1", "Sc1ccccc1", "[S-]c1ccccc1"),
        ("CC(=O)[S-]", "CC(=O)S", "CC(=O)[S-]"),
    ],
    ids=["thiol", "phenyl_thiol", "thioic_acid"],
)
def test_thiolate_input_is_neutralized(
    smiles: str, protonated_form: str, deprotonated_form: str
) -> None:
    """Checks that sulfur given already deprotonated is still recognized as a
    site. The thiol patterns require S-H, so [S-] inputs kept their charge at
    any pH."""

    assert protonated(smiles, VERY_ACIDIC_PH) == [canonical(protonated_form)]
    assert protonated(smiles, VERY_BASIC_PH) == [canonical(deprotonated_form)]


CHARGE_ON_NON_H_NITROGEN = [
    # [charged input, neutral form, id]
    ["C[n+]1cc[nH]c1", "Cn1ccnc1", "imidazolium"],
    ["C[n+]1ccc[nH]1", "Cn1cccn1", "pyrazolium"],
    ["CC(N)=[N+](C)C", "CC(=N)N(C)C", "amidinium"],
    ["CN(C)C(N)=[N+](C)C", "CN(C)C(=N)N(C)C", "guanidinium"],
]


@pytest.mark.parametrize(
    ("charged", "neutral"),
    [group[:2] for group in CHARGE_ON_NON_H_NITROGEN],
    ids=[group[2] for group in CHARGE_ON_NON_H_NITROGEN],
)
def test_charge_on_non_h_nitrogen_is_neutralized(charged: str, neutral: str) -> None:
    """Checks that a cation drawn with the charge on a nitrogen that has no H
    is neutralized. Only N+-H lost its proton, so these resonance forms kept
    their charge and were never matched by the neutral site patterns."""

    assert neutralized(charged) == canonical(neutral)


@pytest.mark.parametrize(
    ("charged", "neutral"),
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
    output."""

    assert protonated(charged, ph) == protonated(neutral, ph)


@pytest.mark.parametrize(
    "smiles",
    ["C[n+]1ccn(C)c1", "CC(N(C)C)=[N+](C)C"],
    ids=["dimethylimidazolium", "peralkyl_amidinium"],
)
def test_fully_substituted_cation_stays_charged(smiles: str) -> None:
    """Checks that the charge-shift rules leave a cation with no N-H alone,
    since it has no proton to lose and really is permanently charged."""

    assert neutralized(smiles) == canonical(smiles)


@pytest.mark.parametrize(
    "smiles",
    [
        "CCOC(=O)C=[N+]=[N-]",
        "[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)O",
        "[N-]=[N+]=CC(=O)CC[C@H](N)C(=O)O",
    ],
    ids=["ethyl_diazoacetate", "azaserine", "DON"],
)
def test_diazo_input_is_not_neutralized(smiles: str) -> None:
    """Checks that the terminal N- of a diazo group is left alone. The rule
    that turns azide N- into N-H also matched it, giving C=[N+]=N, which no
    site pattern deprotonates."""

    assert neutralized(smiles) == canonical(smiles)


@pytest.mark.parametrize(
    "ph", [VERY_ACIDIC_PH, 7.4, VERY_BASIC_PH], ids=["very_acidic", "7.4", "very_basic"]
)
def test_diazo_without_sites_is_unchanged(ph: float) -> None:
    """Checks that a diazo compound with no ionizable site comes out neutral
    at every pH, rather than as a permanent +1 cation."""

    smiles = "CCOC(=O)C=[N+]=[N-]"
    assert protonated(smiles, ph) == [canonical(smiles)]


@pytest.mark.parametrize(
    ("ph", "expected"),
    [
        (VERY_ACIDIC_PH, "[N-]=[N+]=CC(=O)OC[C@H]([NH3+])C(=O)O"),
        (VERY_BASIC_PH, "[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)[O-]"),
    ],
    ids=["very_acidic", "very_basic"],
)
def test_diazo_does_not_shift_net_charge_of_other_sites(
    ph: float, expected: str
) -> None:
    """Checks that azaserine's amine and acid are protonated as usual while
    its diazo group stays neutral, so the net charge is +1 and -1 rather than
    +2 and 0."""

    assert protonated("[N-]=[N+]=CC(=O)OC[C@H](N)C(=O)O", ph) == [canonical(expected)]


def test_azide_input_is_still_neutralized() -> None:
    """Checks that the diazo exception does not reach azides, whose terminal
    N- must still gain an H for the Azide site to match."""

    assert neutralized("CN=[N+]=[N-]") == canonical("CN=[N+]=N")


def test_neutralization_raises_when_a_rule_never_converges() -> None:
    """Checks that a rule whose product still matches its reactant raises,
    naming the rule and the input. The loop had no pass limit, so one such
    rule hung the whole batch."""

    neutralizer = MoleculeNeutralizer(rxn_data=(("[Sv1-1:1]", "[S-1:1]"),))

    with pytest.raises(RuntimeError, match="did not converge") as info:
        neutralizer.neutralize_smiles("CC[S-]")
    assert "[Sv1-1:1] >> [S-1:1]" in str(info.value)


def test_neutralization_enumerates_one_product() -> None:
    """Checks that each pass asks RDKit for a single product. Only the first
    is used, but the default enumerated up to 1000 per pass."""

    requested: list[int] = []

    class RecordingReaction:
        """Wraps a reaction to record how many products each call asks for,
        since RDKit's C++ methods cannot be patched directly."""

        def __init__(self, rxn: rdChemReactions.ChemicalReaction) -> None:
            """Stores the real reaction.

            Args:
                rxn: The compiled reaction to delegate to.
            """

            self.rxn = rxn

        def RunReactants(
            self, reactants: tuple[Chem.Mol, ...], maxProducts: int = 1000
        ) -> tuple[tuple[Chem.Mol, ...], ...]:
            """Records maxProducts, then runs the real reaction.

            Args:
                reactants: Passed through unchanged.
                maxProducts: The product limit to record.

            Returns:
                The real reaction's products.
            """

            requested.append(maxProducts)
            return self.rxn.RunReactants(reactants, maxProducts)

    neutralizer = MoleculeNeutralizer()
    for reaction in neutralizer.registry.reactions:
        reaction._rxn = RecordingReaction(reaction._rxn)

    smiles = neutralizer.neutralize_smiles("[O-]C(=O)CC(=O)[O-]")

    assert canonical(smiles) == canonical("OC(=O)CC(=O)O")
    assert requested != [] and all(n == 1 for n in requested), requested
