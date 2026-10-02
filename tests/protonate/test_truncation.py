"""Checks how max_variants chooses which protonation states to keep."""

from collections import Counter

import pytest
from loguru import logger
from rdkit import Chem

import dimorphite_dl.protonate.run as run_module
from dimorphite_dl import protonate_smiles
from dimorphite_dl.protonate.change import charge_log_probabilities
from dimorphite_dl.protonate.site import ProtonationSite, ProtonationState

# Wide enough that every carboxyl and amine site is BOTH at pH 7.
EVERY_SITE_BOTH_PRECISION = 10.0

# Eight carboxyl groups on distinct, non-equivalent carbons, so every one of
# the 2**8 charge patterns is a different molecule.
EIGHT_CARBOXYLS = "CC(C(=O)O)" * 8
EIGHT_CARBOXYL_STATES = 2**8


def canonical(smiles: str) -> str:
    """Canonicalizes so expected values can be written in any valid form.

    Args:
        smiles: A valid SMILES string.

    Returns:
        RDKit's canonical isomeric SMILES for smiles.
    """

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


def formal_charge(smiles: str) -> int:
    """Sums a molecule's formal charges, which for the carboxyl tests is minus
    the number of deprotonated groups.

    Args:
        smiles: A valid SMILES string.

    Returns:
        The molecule's net formal charge.
    """

    return sum(atom.GetFormalCharge() for atom in Chem.MolFromSmiles(smiles).GetAtoms())


def protonate_at_7(smiles: str, max_variants: int) -> list[str]:
    """Runs the pipeline at pH 7 with every site BOTH.

    Args:
        smiles: The input SMILES string.
        max_variants: The most states to return.

    Returns:
        The output SMILES.
    """

    return protonate_smiles(
        smiles,
        ph_min=7.0,
        ph_max=7.0,
        precision=EVERY_SITE_BOTH_PRECISION,
        max_variants=max_variants,
    )


@pytest.mark.parametrize(
    ("pka", "ph"), [(7.0, 7.0), (4.0, 7.0), (10.0, 7.0), (-1000.0, 7.4)]
)
def test_charge_log_probabilities_are_fractions(pka: float, ph: float) -> None:
    """Checks that a BOTH site's deprotonated and protonated fractions sum to
    one, split evenly at the pKa, and do not overflow at Nitro's pKa."""

    deprotonated, protonated = charge_log_probabilities(ProtonationState.BOTH, pka, ph)

    assert 10**deprotonated + 10**protonated == pytest.approx(1.0)
    if pka == ph:
        assert deprotonated == pytest.approx(protonated)
    else:
        assert (protonated > deprotonated) == (pka > ph)


@pytest.mark.parametrize(
    "state", [ProtonationState.PROTONATED, ProtonationState.DEPROTONATED]
)
def test_single_state_site_does_not_affect_ranking(state: ProtonationState) -> None:
    """Checks that a site with only one charge scores zero, since it shifts
    every state equally."""

    assert charge_log_probabilities(state, 8.0, 7.0) == [0.0]


@pytest.mark.parametrize(
    ("max_variants", "expected"),
    [
        (1, ["[NH3+]CC(=O)[O-]"]),
        (2, ["[NH3+]CC(=O)[O-]", "NCC(=O)[O-]"]),
        (3, ["[NH3+]CC(=O)[O-]", "NCC(=O)[O-]", "[NH3+]CC(=O)O"]),
    ],
    ids=["one", "two", "three"],
)
def test_max_variants_keeps_most_probable_states(
    max_variants: int, expected: list[str]
) -> None:
    """Checks that truncation keeps the likeliest states at the pH. At pH 7
    the amine (pKa about 8.2) is mostly protonated and the carboxyl (about
    3.5) mostly deprotonated. Taking the first states generated instead kept
    whatever the enumeration happened to produce first.

    Glycine rather than a longer amino acid, because its amine sits next to
    the carboxyl and so falls to Amines_primary_secondary_tertiary, whose mean
    of 8.16 is what the ranking above assumes. An amine with a plain alkyl
    chain is claimed by Alkylamine_primary instead, and its mean of 10.55
    leaves the two sites' states nearly equally likely, which would make the
    expected order here a coin flip."""

    output = protonate_at_7("NCC(=O)O", max_variants)

    assert sorted(canonical(s) for s in output) == sorted(
        canonical(s) for s in expected
    ), output


def test_max_variants_does_not_pin_later_sites() -> None:
    """Checks that with eight equivalent carboxyls at pH 7 and room for nine
    states, the cap keeps the fully deprotonated state and all eight singly
    protonated ones. Enumeration stopped as soon as the cap was reached, so
    every later carboxyl kept its input state."""

    output = protonate_at_7(EIGHT_CARBOXYLS, 9)

    assert sorted(formal_charge(s) for s in output) == [-8] + [-7] * 8, output


def test_many_both_sites_uncapped() -> None:
    """Checks that every combination is produced when the cap is not hit, and
    that no truncation warning is logged."""

    messages: list[str] = []
    handler = logger.add(lambda message: messages.append(str(message)), level="WARNING")
    try:
        output = protonate_at_7(EIGHT_CARBOXYLS, EIGHT_CARBOXYL_STATES)
    finally:
        logger.remove(handler)

    assert len(output) == EIGHT_CARBOXYL_STATES
    assert len(set(output)) == EIGHT_CARBOXYL_STATES
    assert all(Chem.MolFromSmiles(s) is not None for s in output)
    assert not any("Limited number of variants" in m for m in messages), messages


@pytest.mark.parametrize("max_variants", [1, 10, 128], ids=["one", "ten", "default"])
def test_many_both_sites_capped(
    max_variants: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Checks that the cap limits the output, that every site is still
    applied after the cap is reached, that no site is applied to more than
    max_variants parents (so memory stays bounded), and that the truncation
    is logged."""

    original = run_module.protonate_site
    parents_per_site: Counter[tuple[int, ...]] = Counter()

    def recording_protonate_site(
        mols: list[Chem.Mol],
        site: ProtonationSite,
        ph_min: float,
        ph_max: float,
        precision: float,
    ) -> list[Chem.Mol]:
        """Counts how many parent states each site is applied to.

        Args:
            mols: The states going into this call.
            site: The site being protonated.
            ph_min: Passed through unchanged.
            ph_max: Passed through unchanged.
            precision: Passed through unchanged.

        Returns:
            The original function's result.
        """

        parents_per_site[site.idxs_match] += len(mols)
        return original(mols, site, ph_min, ph_max, precision)

    monkeypatch.setattr(run_module, "protonate_site", recording_protonate_site)

    messages: list[str] = []
    handler = logger.add(lambda message: messages.append(str(message)), level="WARNING")
    try:
        output = protonate_at_7(EIGHT_CARBOXYLS, max_variants)
    finally:
        logger.remove(handler)

    assert len(output) == max_variants
    assert len(set(output)) == max_variants
    assert all(Chem.MolFromSmiles(s) is not None for s in output)

    assert len(parents_per_site) == 8
    assert max(parents_per_site.values()) <= max_variants

    assert any("Limited number of variants" in m for m in messages), messages
