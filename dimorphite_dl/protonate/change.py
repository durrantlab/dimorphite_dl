import copy

from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol

from dimorphite_dl.protonate.site import ProtonationSite, ProtonationState


def protonate_site(
    mols: list[Mol], site: ProtonationSite, ph_min, ph_max, precision
) -> list[Mol]:
    """Protonate a specific site in a list of molecules.

    Args:
        mols: List of molecule objects
        site: ProtonationSite object with protonation information

    Returns:
        List of appropriately protonated molecule objects
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
            current_mols = branched
    return current_mols


def set_protonation_charge(
    mols: list[Mol], idx: int, charges: list[int], prot_site_name: str
) -> list[Mol]:
    """Set atomic charge on a specific site for a set of molecules.

    Args:
        mols: List of input molecule objects
        idx: Index of the atom to modify
        charges: List of charges to assign at this site
        prot_site_name: Name of the protonation site

    Returns:
        List of processed molecule objects
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
            except Exception as e:
                logger.warning(
                    "Error processing molecule with charge {}: {}", charge, str(e)
                )
                continue
            mols_charged.append(processed_mol)
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

        # Special case for aromatic nitrogen deprotonation
        mol_smiles = Chem.MolToSmiles(mol_copy)
        if "[nH-]" in mol_smiles:
            atom.SetNumExplicitHs(0)

        # Update property cache
        mol_copy.UpdatePropertyCache(strict=False)

    except Exception as e:
        logger.warning("Error setting atom properties: {}", str(e))
        return None

    return mol_copy


def _set_nitrogen_properties(
    atom: Chem.Atom, charge: int, bond_order_total: int
) -> None:
    """Set properties for nitrogen atoms based on charge and bonding."""
    atom_idx = atom.GetIdx()
    is_aromatic = atom.GetIsAromatic()
    degree = atom.GetDegree()
    logger.debug(
        "Setting N properties: index={}, charge={}, bond_order={}, aromatic={}, degree={}",
        atom_idx,
        charge,
        bond_order_total,
        is_aromatic,
        degree,
    )

    # Don't protonate aromatic nitrogens that already have 3 neighbors
    # This handles cases like the N9 in adenine that's connected to a sugar
    if is_aromatic and degree == 3 and charge > 0:
        logger.debug("N is aromatic, with a degree of 3, cannot set positive charge")
        # This nitrogen is already saturated, don't change it
        return
    atom.SetFormalCharge(charge)

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

    h_count = h_count_map.get((charge, bond_order_total), -1)
    if h_count != -1:
        atom.SetNumExplicitHs(h_count)


def _set_other_element_properties(
    atom: Chem.Atom, charge: int, element: int, bond_order_total: float
) -> None:
    """Set properties for non-nitrogen atoms."""
    atom.SetFormalCharge(charge)

    # Special handling for oxygen and sulfur
    if element in (8, 16):  # O and S
        if charge == 0 and bond_order_total == 1:
            atom.SetNumExplicitHs(1)
        elif charge == -1 and bond_order_total == 1:
            atom.SetNumExplicitHs(0)
