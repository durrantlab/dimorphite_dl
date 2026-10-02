import copy
import math

from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol

from dimorphite_dl.protonate.site import ProtonationSite, ProtonationState


def protonate_site(
    mols: list[Mol],
    site: ProtonationSite,
    ph_min: float,
    ph_max: float,
    precision: float,
) -> list[Mol]:
    """Protonate a specific site in a list of molecules.

    Args:
        mols: List of molecule objects.
        site: ProtonationSite object with protonation information.
        ph_min: Minimum pH to expose the site to.
        ph_max: Maximum pH to expose the site to.
        precision: pKa standard deviation prefactor to consider.

    Returns:
        List of appropriately protonated molecule objects. If there is any issue,
            this will return an empty list.
    """
    if not mols:
        logger.warning("No molecules provided for protonation")
        return []

    logger.debug("Protonating site: {}", site.name)

    unique_states = list(site.get_unique_states(ph_min, ph_max, precision))

    current_mols = mols

    for idx_atom, state in unique_states:
        charges = state.get_charges()

        # If the state is not BOTH, we apply its single charge to each
        # molecule in current_mols without creating branches.
        if state != ProtonationState.BOTH:
            logger.debug(
                "Site {} atom {} has exclusive state {}; applying to all molecules",
                site.name,
                idx_atom,
                state.to_str(),
            )
            processed = set_protonation_charge(
                current_mols, idx_atom, charges, site.name
            )
            if len(processed) == 0:
                return []
            current_mols = processed

        else:
            logger.debug(
                "Site {} atom {} is BOTH; branching into {} variants per molecule",
                site.name,
                idx_atom,
                charges,
            )

            branched = []
            for mol in current_mols:
                try:
                    variants = set_protonation_charge(
                        [mol], idx_atom, charges, site.name
                    )
                    branched.extend(variants)
                except Exception as e:
                    logger.error("Error protonating site {}: {}", idx_atom, str(e))
                    return []
            current_mols = branched
    return current_mols


def log10_one_plus_pow10(x: float) -> float:
    """Computes log10(1 + 10**x) without evaluating 10**x, which overflows
    for the extreme pKa values in the SMARTS file (e.g., Nitro's -1000).

    Args:
        x: The exponent.

    Returns:
        log10(1 + 10**x).
    """
    return max(x, 0.0) + math.log10(1.0 + 10.0 ** -abs(x))


def charge_log_probabilities(
    state: ProtonationState, pka: float, ph: float
) -> list[float]:
    """Scores each charge a pKa site is enumerated with, so that truncation
    can keep the likeliest states rather than the first ones generated.

    Uses the Henderson-Hasselbalch fraction at the site's mean pKa. A state
    with only one charge shifts every variant's score equally, so it scores
    0 rather than a large constant that would cost precision.

    Args:
        state: The site's target state at the pH range.
        pka: The site's mean pKa.
        ph: The pH to score at.

    Returns:
        log10 of the fraction in each charge state, in the order of
        state.get_charges().
    """
    charges = state.get_charges()
    if len(charges) == 1:
        return [0.0]

    # Protonated fraction: 1 / (1 + 10**(pH - pKa)). Deprotonated fraction:
    # 1 / (1 + 10**(pKa - pH)).
    return [
        (
            -log10_one_plus_pow10(ph - pka)
            if charge == 0
            else -log10_one_plus_pow10(pka - ph)
        )
        for charge in charges
    ]


def site_log_probabilities(
    site: ProtonationSite, ph_min: float, ph_max: float, precision: float, ph: float
) -> list[float]:
    """Scores every variant protonate_site makes from one parent molecule,
    in the same order, so callers can rank the variants it returns.

    Args:
        site: The protonation site.
        ph_min: Minimum pH of the range, which sets each pKa's state.
        ph_max: Maximum pH of the range.
        precision: pKa standard deviation prefactor.
        ph: The pH to score at.

    Returns:
        log10 probability of each variant, relative to its parent.
    """
    # Mirrors get_unique_states: one entry per distinct (atom, state), in
    # pKa order, keeping the first pKa that produced it.
    pka_by_state: dict[tuple[int, ProtonationState], float] = {}
    for pka in site.pkas:
        key = (site.idxs_match[pka.idx_site], pka.get_state(ph_min, ph_max, precision))
        pka_by_state.setdefault(key, pka.mean)

    # protonate_site branches parent-major: each existing variant is followed
    # by its copies for each charge in turn.
    scores = [0.0]
    for (_, state), mean in pka_by_state.items():
        charge_scores = charge_log_probabilities(state, mean, ph)
        scores = [parent + charge for parent in scores for charge in charge_scores]
    return scores


def set_protonation_charge(
    mols: list[Mol], idx: int, charges: list[int], prot_site_name: str
) -> list[Mol]:
    """Set atomic charge on a specific site for a set of molecules.

    Args:
        mols: List of input molecule objects.
        idx: Index of the atom to modify.
        charges: List of charges to assign at this site.
        prot_site_name: Name of the protonation site.

    Returns:
        List of processed molecule objects. If anything goes wrong, then we return
            an empty list.
    """
    is_special_nitrogen = "*" in prot_site_name

    mols_charged = []
    for charge in charges:
        nitrogen_charge = charge + 1

        # Special case for nitrogen moieties where acidic group is neutral
        if is_special_nitrogen:
            nitrogen_charge = nitrogen_charge - 1

        for mol in mols:
            try:
                processed_mol = _apply_charge_to_molecule(
                    mol, idx, charge, nitrogen_charge
                )
                if processed_mol is not None:
                    mols_charged.append(processed_mol)
                else:
                    return []
            except Exception as e:
                logger.warning(
                    "Error processing molecule with charge {}: {}", charge, str(e)
                )
                return []
    return mols_charged


def _apply_charge_to_molecule(
    mol: Mol, idx: int, charge: int, nitrogen_charge: int
) -> Mol | None:
    """Apply charge to a specific atom in a molecule.

    Args:
        mol: Input molecule
        idx: Atom index
        charge: Charge for non-nitrogen atoms
        nitrogen_charge: Charge for nitrogen atoms
        prot_site_name: Name of protonation site

    Returns:
        Modified molecule or None if processing fails
    """
    logger.trace(
        "Applying charge of {} at index {} to SMILES: {}",
        charge,
        idx,
        Chem.MolToSmiles(mol),
    )
    # Create deep copy to avoid modifying original
    mol_copy = copy.deepcopy(mol)

    # Remove hydrogens first
    try:
        mol_copy = Chem.RemoveHs(mol_copy)
        if mol_copy is None:
            logger.warning("RemoveHs returned None for molecule")
            return None
    except Exception as e:
        logger.warning("Failed to remove hydrogens: {}", str(e))
        return None

    # Kept so an unbuildable state can fall back to the parent's state.
    mol_parent = Chem.Mol(mol_copy)

    # Validate atom index
    if idx >= mol_copy.GetNumAtoms():
        logger.warning(
            "Atom index {} out of range (molecule has {} atoms)",
            idx,
            mol_copy.GetNumAtoms(),
        )
        return None

    atom = mol_copy.GetAtomWithIdx(idx)
    element = atom.GetAtomicNum()

    # Calculate explicit bond order
    try:
        explicit_bond_order_total = sum(
            b.GetBondTypeAsDouble() for b in atom.GetBonds()
        )
    except Exception as e:
        logger.warning("Error calculating bond order for atom {}: {}", idx, str(e))
        return None

    # Set formal charge and explicit hydrogens based on element type
    try:
        if element == 7:  # Nitrogen
            _set_nitrogen_properties(atom, nitrogen_charge, explicit_bond_order_total)
        else:
            _set_other_element_properties(
                atom, charge, element, explicit_bond_order_total
            )

        # A deprotonated aromatic N-H keeps its H here, because the bond-order
        # table has no entry for two aromatic bonds. This was detected by
        # searching the SMILES for "[nH-]", which misses atoms written with a
        # map number or isotope ("[nH-:1]", "[15nH-]"), so those sites could
        # never be deprotonated.
        if (
            element == 7
            and atom.GetIsAromatic()
            and atom.GetFormalCharge() == -1
            and atom.GetNumExplicitHs() > 0
        ):
            logger.debug("Aromatic N- still carries an H; setting its H count to 0")
            atom.SetNumExplicitHs(0)

        # Update property cache
        mol_copy.UpdatePropertyCache(strict=False)

    except Exception as e:
        logger.warning("Error setting atom properties: {}", str(e))
        return None

    # Some states have no valid structure when only this atom is edited (e.g.,
    # a bridgehead aromatic N given a +1 charge cannot be kekulized). The site
    # keeps its parent's state; otherwise the invalid SMILES is rejected at
    # the end, and when it is the site's only state the molecule disappears.
    sanitized = Chem.SanitizeMol(Chem.Mol(mol_copy), catchErrors=True)
    if sanitized != Chem.SanitizeFlags.SANITIZE_NONE:
        logger.warning(
            "No valid state with charge {} at atom {} of {}; keeping the "
            "unmodified structure in its place",
            charge,
            idx,
            Chem.MolToSmiles(mol_parent),
        )
        return mol_parent

    return mol_copy


def _set_nitrogen_properties(
    atom: Chem.Atom, charge: int, bond_order_total: int
) -> None:
    """Set properties for nitrogen atoms based on charge and bonding."""
    atom_idx = atom.GetIdx()
    is_aromatic = atom.GetIsAromatic()
    degree = atom.GetDegree()
    logger.trace(
        "Setting N properties: index={}, charge={}, bond_order={}, aromatic={}, degree={}",
        atom_idx,
        charge,
        bond_order_total,
        is_aromatic,
        degree,
    )

    # Handling niche cases of aromatics often detected on NADP
    if charge == 1 and bond_order_total == 4.0 and is_aromatic and degree == 3:
        return

    atom.SetFormalCharge(charge)
    logger.debug("Set formal charge to {}", charge)

    # Set explicit hydrogens based on charge and bond order
    h_count_map = {
        (1, 1): 3,
        (1, 2): 2,
        (1, 3): 1,  # Positive charge
        (0, 1): 2,
        (0, 2): 1,  # Neutral
        (-1, 1): 1,
        (-1, 2): 0,  # Negative charge
    }

    h_count = h_count_map.get((int(charge), int(bond_order_total)), -1)
    if h_count != -1:
        logger.debug("Setting hydrogen count to {}", h_count)
        atom.SetNumExplicitHs(h_count)


def _set_other_element_properties(
    atom: Chem.Atom, charge: int, element: int, bond_order_total: float
) -> None:
    """Set properties for non-nitrogen atoms."""
    atom_idx = atom.GetIdx()
    is_aromatic = atom.GetIsAromatic()
    degree = atom.GetDegree()
    logger.trace(
        "Setting {} properties: index={}, charge={}, bond_order={}, aromatic={}, degree={}",
        element,
        atom_idx,
        charge,
        bond_order_total,
        is_aromatic,
        degree,
    )

    atom.SetFormalCharge(charge)
    logger.debug("Set formal charge to {}", charge)

    # Special handling for oxygen and sulfur
    if element in (8, 16):  # O and S
        if charge == 0 and bond_order_total == 1:
            atom.SetNumExplicitHs(1)
            logger.debug("Set explicit hydrogens for this atom to 1")
        elif charge == -1 and bond_order_total == 1:
            atom.SetNumExplicitHs(0)
            logger.debug("Set explicit hydrogens for this atom to 0")