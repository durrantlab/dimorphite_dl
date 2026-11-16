"""Checks the site entries whose pKas come from the literature rather than the
training data: Hydrazoic_acid, Carboxamide, Hydrazide, Acylhydrazone,
Diazine, and Tetrazole_NH."""

import pytest
from rdkit import Chem

from dimorphite_dl import protonate_smiles
from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.protonate.detect import ProtonationSiteDetector

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


def claiming_entries(smiles: str) -> set[str]:
    """Reports which entries claim a site, so that the scope of a new entry can
    be tested without asserting an end state for compounds whose pKa the entry
    does not stand behind.

    Args:
        smiles: The input SMILES string.

    Returns:
        The names of the site entries that matched, as they appear in
        site_substructures.smarts.
    """
    mol = MoleculeRecord(smiles).prepare_for_protonation()
    detector = ProtonationSiteDetector()
    return {site.name for site in detector._detect_all_sites_in_molecule(mol)}


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


def test_isoniazid_hydrazide_is_neutral_at_physiological_ph() -> None:
    """The headline case for Hydrazide: Amide_electronegative assigned the
    hydrazide N-H a pKa of 3.49, so isoniazid was returned only as the anion.
    The measured value is 10.95, and the pyridine is not a site at pH 7.4
    either, so the neutral molecule is the only state here."""
    output = protonate_at("NNC(=O)c1ccncc1", PHYSIOLOGICAL_PH)
    assert canonical_set(output) == canonical_set(["NNC(=O)c1ccncc1"]), output


def test_isoniazid_hydrazide_is_neutral_at_default_ph_range() -> None:
    """Isoniazid also gives a pyridinium variant over the default 6.4 to 8.4
    window, because the Aromatic_nitrogen_unprotonated window reaches 6.425
    and so overlaps ph_min by 0.025. That is the generic aromatic nitrogen
    rule, not Hydrazide, and it predates this entry, so only the hydrazide
    nitrogen is asserted here."""
    for out in protonate_smiles("NNC(=O)c1ccncc1"):
        mol = Chem.MolFromSmiles(out)
        charged_n = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "N" and atom.GetFormalCharge() < 0
        ]
        assert not charged_n, out


def test_acetohydrazide_is_neutral_at_default_ph_range() -> None:
    """The same claim on a hydrazide with no other site at all, so the whole
    output can be pinned. The default window is what the command line runs,
    and it is wider than any single pH, so a mean this entry got wrong would
    show up here as an extra anionic variant."""
    output = protonate_smiles("CC(=O)NN")
    assert canonical_set(output) == canonical_set(["CC(=O)NN"]), output


def test_trifluoroacetohydrazide_is_anion_at_default_ph_range() -> None:
    """The companion check to the one above, and the reason Hydrazide carries
    a halogen screen: the perfluorinated hydrazides are measured at 7.92 to
    8.02, so the anion Amide_electronegative gives them is roughly right and
    should survive at the default window."""
    output = protonate_smiles("NNC(=O)C(F)(F)F")
    assert canonical_set(output) == canonical_set(["N[N-]C(=O)C(F)(F)F"]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=O)NN", id="acetohydrazide_13.04"),
        pytest.param("COCC(=O)NN", id="methoxyacetohydrazide_12.14"),
        pytest.param("N#CCC(=O)NN", id="cyanoacetohydrazide_11.04"),
        pytest.param("NNC(=O)c1ccncc1", id="isoniazid_10.8"),
        pytest.param("NNC(=O)c1ccccc1", id="benzohydrazide"),
        pytest.param("CNNC(C)=O", id="N-methylacetohydrazide"),
        pytest.param("O=C(NN)c1ccc(Cl)cc1", id="4-chlorobenzohydrazide"),
    ],
)
def test_hydrazide_nh_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Hydrazide N-H runs from about 11 to 13, so no hydrazide anion should
    appear in any variant. Other sites (e.g., the second nitrogen as a base)
    may still vary, so only N- is checked.

    Args:
        smiles: A molecule with a hydrazide N-H in scope for Hydrazide.
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
        pytest.param("CC(=O)NN", id="acetohydrazide"),
        pytest.param("NNC(=O)c1ccncc1", id="isoniazid"),
        pytest.param("CNNC(C)=O", id="N-methylacetohydrazide"),
    ],
)
def test_hydrazide_claims_the_site_from_amide_electronegative(smiles: str) -> None:
    """The fix is a routing change, so check the routing directly: these N-H
    were claimed by Amide_electronegative and should now be claimed by
    Hydrazide instead.

    Args:
        smiles: A molecule with a hydrazide N-H in scope for Hydrazide.
    """
    entries = claiming_entries(smiles)
    assert "*Hydrazide" in entries, entries
    assert "*Amide_electronegative" not in entries, entries


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=C(NN=Cc1ccccc1)c1ccccc1", id="benzaldehyde_benzoylhydrazone"),
        pytest.param("CC(=O)NN=Cc1ccccc1", id="benzaldehyde_acetohydrazone"),
        pytest.param("O=C(NN=C(C)C)c1ccccc1", id="acetone_benzoylhydrazone"),
        pytest.param(
            "O=C(NN=Cc1ccccc1)c1ccncc1", id="isoniazid_benzaldehyde_hydrazone"
        ),
    ],
)
def test_acylhydrazone_nh_stays_neutral_at_physiological_ph(smiles: str) -> None:
    """Acylhydrazone N-H measures near 11.6, so no hydrazone anion should
    appear in any variant. Amide_electronegative gave this class 3.49 and so
    returned the anion alone. Other sites may still vary, so only N- is
    checked.

    Args:
        smiles: A molecule with an acylhydrazone N-H in scope.
    """
    for out in protonate_at(smiles, PHYSIOLOGICAL_PH):
        mol = Chem.MolFromSmiles(out)
        charged_n = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "N" and atom.GetFormalCharge() < 0
        ]
        assert not charged_n, out


def test_benzoylhydrazone_is_neutral_at_default_ph_range() -> None:
    """The same claim on an acylhydrazone with no other site, so the whole
    output can be pinned over the default 6.4 to 8.4 window."""
    output = protonate_smiles("O=C(NN=Cc1ccccc1)c1ccccc1")
    assert canonical_set(output) == canonical_set(["O=C(NN=Cc1ccccc1)c1ccccc1"]), output


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("O=C(NN=Cc1ccccc1)c1ccccc1", id="benzoylhydrazone"),
        pytest.param("CC(=O)NN=Cc1ccccc1", id="acetohydrazone"),
    ],
)
def test_acylhydrazone_claims_the_site_from_amide_electronegative(smiles: str) -> None:
    """The routing half of the fix, as for Hydrazide: these N-H were claimed
    by Amide_electronegative and should now be claimed by Acylhydrazone.

    Args:
        smiles: A molecule with an acylhydrazone N-H in scope.
    """
    entries = claiming_entries(smiles)
    assert "*Acylhydrazone" in entries, entries
    assert "*Amide_electronegative" not in entries, entries


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("CC(=O)NN", id="acetohydrazide"),
        pytest.param("NNC(=O)c1ccncc1", id="isoniazid"),
    ],
)
def test_hydrazide_and_acylhydrazone_do_not_overlap(smiles: str) -> None:
    """The two entries are mutually exclusive by hybridization of the second
    nitrogen, NX3 against NX2, so a plain hydrazide must not reach the
    hydrazone entry whatever order they sit in.

    Args:
        smiles: A plain hydrazide.
    """
    entries = claiming_entries(smiles)
    assert "*Hydrazide" in entries, entries
    assert "*Acylhydrazone" not in entries, entries


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NNC(=O)C(F)(F)F", id="trifluoroacetohydrazide"),
        pytest.param("NNC(=O)C(F)F", id="difluoroacetohydrazide"),
        pytest.param("NC(=O)NN", id="semicarbazide"),
        pytest.param("CC(=O)NNC(C)=O", id="1,2-diacetylhydrazine"),
        pytest.param("CCOC(=O)NN", id="ethyl_carbazate"),
        pytest.param("NNS(=O)(=O)c1ccccc1", id="benzenesulfonohydrazide"),
        pytest.param("FC(F)(F)C(=O)NN=Cc1ccccc1", id="trifluoroacetyl_hydrazone"),
        pytest.param("NC(=O)NN=Cc1ccccc1", id="semicarbazone"),
        pytest.param("NC(=S)NN=Cc1ccccc1", id="thiosemicarbazone"),
        pytest.param("O=C(NN=C(N)c1ccccc1)c1ccccc1", id="amidrazone"),
    ],
)
def test_classes_excluded_from_the_new_entries_keep_their_entry(smiles: str) -> None:
    """Both new entries are scoped to what their measurements cover. The
    perfluorinated hydrazides sit at 7.92 to 9.66, and no aqueous series was
    found for the remaining classes, so none of them should be claimed by
    either entry. Routing is asserted rather than an end state, since the pKa
    these keep is not one either entry stands behind.

    Args:
        smiles: A molecule neither Hydrazide nor Acylhydrazone claims.
    """
    entries = claiming_entries(smiles)
    assert "*Hydrazide" not in entries, entries
    assert "*Acylhydrazone" not in entries, entries


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


@pytest.mark.parametrize(
    ("smiles", "anion"),
    [
        pytest.param("c1nn[nH]n1", "c1nn[n-]n1", id="tetrazole_4.86"),
        pytest.param("Cc1nn[nH]n1", "Cc1nn[n-]n1", id="5-methyltetrazole_5.50"),
        pytest.param("c1ccc(-c2nn[nH]n2)cc1", "c1ccc(-c2nn[n-]n2)cc1", id="5-phenyl"),
    ],
)
def test_tetrazole_is_anion_at_physiological_ph(smiles: str, anion: str) -> None:
    """Tetrazole N-H has a pKa near 5, so the anion is the only state at pH
    7.4. The generic aromatic N-H rule returned both.

    Args:
        smiles: A neutral 5-substituted tetrazole.
        anion: Its N-deprotonated form.
    """
    output = protonate_at(smiles, PHYSIOLOGICAL_PH)
    assert canonical_set(output) == canonical_set([anion]), output


def test_tetrazole_not_amide_is_deprotonated() -> None:
    """The case from a GitHub issue: with max_variants=1 an amide anion was
    returned and the tetrazole left neutral. The amide should stay neutral and the
    tetrazole should be the only site that ionizes, in every variant."""
    output = protonate_smiles(
        "O=C(NCc1cccc(-c2nnn[nH]2)c1)c1ccccc1",
        ph_min=6.9,
        ph_max=7.5,
        precision=1.0,
    )
    assert canonical_set(output) == canonical_set(
        ["O=C(NCc1cccc(-c2nnn[n-]2)c1)c1ccccc1"]
    ), output


def test_aminotetrazole_keeps_generic_rule() -> None:
    """5-Aminotetrazole (pKa about 6.1) is excluded from Tetrazole_NH, so it
    still gets the generic window and comes back in both states."""
    output = protonate_at("Nc1nn[nH]n1", PHYSIOLOGICAL_PH)
    assert len(output) > 1, output
