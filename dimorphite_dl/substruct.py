from typing import Any

import copy
import importlib.resources as pkg_resources
from dataclasses import dataclass

from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol

from dimorphite_dl import data  # type: ignore
from dimorphite_dl import protect
from dimorphite_dl.mol import smiles_to_mol


@dataclass
class ProtonationSite:
    """Data class for protonation site information."""

    atom_idx: int
    target_state: str
    site_name: str

    def __post_init__(self):
        if self.target_state not in {"PROTONATED", "DEPROTONATED", "BOTH"}:
            raise ValueError(f"Invalid target state: {self.target_state}")


@dataclass
class SubstructureData:
    """Data class for substructure information."""

    name: str
    smart: str
    mol: Mol | None
    prot_states_for_pH: list[list[str]]

    def __post_init__(self):
        if self.mol is None:
            logger.warning(
                "Failed to create mol object for SMARTS: {} ({})", self.smart, self.name
            )


def load_substructure_smarts_file() -> list[str]:
    """Load the substructure SMARTS file, filtering out comments and blank lines.

    Returns:
        List of valid SMARTS lines from the file.

    Raises:
        FileNotFoundError: If the substructure file cannot be found.
        IOError: If there are issues reading the file.
    """
    logger.trace("Loading substructure data from site_substructures.smarts")

    try:
        with pkg_resources.open_text(data, "site_substructures.smarts") as f:
            lines = []
            line_count = 0
            valid_count = 0

            for line in f:
                line_count += 1
                stripped = line.strip()

                # Skip empty lines and comments
                if stripped and not stripped.startswith("#"):
                    lines.append(stripped)
                    valid_count += 1

            logger.debug(
                "Loaded {} valid SMARTS patterns from {} total lines",
                valid_count,
                line_count,
            )
            return lines

    except FileNotFoundError:
        logger.error("Could not find site_substructures.smarts file")
        raise
    except Exception as e:
        logger.error("Error reading substructure file: {}", str(e))
        raise IOError(f"Failed to read substructure file: {e}")


def load_protonation_substructs_calc_state_for_ph(
    min_ph: float = 6.4, max_ph: float = 8.4, pka_std_range: float = 1.0
) -> list[SubstructureData]:
    """Load protonation substructures with calculated states for pH range.

    Args:
        min_ph: Lower bound of pH range (default: 6.4)
        max_ph: Upper bound of pH range (default: 8.4)
        pka_std_range: Precision factor for pKa calculations (default: 1.0)

    Returns:
        List of SubstructureData objects with protonation information.

    Raises:
        ValueError: If pH range is invalid.
    """
    if min_ph > max_ph:
        raise ValueError(f"Invalid pH range: min_ph ({min_ph}) >= max_ph ({max_ph})")

    if pka_std_range < 0:
        raise ValueError(f"pka_std_range must be positive, got: {pka_std_range}")

    logger.info(
        "Loading protonation substructures for pH range {:.1f}-{:.1f}", min_ph, max_ph
    )

    subs = []
    failed_count = 0

    for line_num, line in enumerate(load_substructure_smarts_file(), 1):
        try:
            sub_data = _parse_substructure_line(line, min_ph, max_ph, pka_std_range)
            if sub_data:
                subs.append(sub_data)
            else:
                failed_count += 1
        except Exception as e:
            logger.warning("Failed to parse line {}: '{}' - {}", line_num, line, str(e))
            failed_count += 1

    logger.info(
        "Successfully loaded {} substructures, {} failed", len(subs), failed_count
    )
    return subs


def _parse_substructure_line(
    line: str, min_ph: float, max_ph: float, pka_std_range: float
) -> SubstructureData | None:
    """Parse a single line from the substructure file.

    Args:
        line: Line from the substructure file
        min_ph: Minimum pH
        max_ph: Maximum pH
        pka_std_range: pKa standard range multiplier

    Returns:
        SubstructureData object or None if parsing fails
    """
    parts = line.split()
    if len(parts) < 3:
        logger.warning("Invalid line format (too few parts): '{}'", line)
        return None

    name = parts[0]
    smart = parts[1]
    logger.trace("Parsing SMARTS of {}", name)
    # Create mol object from SMARTS
    try:
        mol = Chem.MolFromSmarts(smart)
        if mol is None:
            logger.warning("Invalid SMARTS pattern: {} ({})", smart, name)
            return None
    except Exception as e:
        logger.warning(
            "Error creating mol from SMARTS '{}' ({}): {}", smart, name, str(e)
        )
        return None

    # Parse pKa ranges (groups of 3: site, mean, std)
    pka_data = parts[2:]
    if len(pka_data) % 3 != 0:
        logger.warning(
            "Invalid pKa data format for {}: expected groups of 3, got {}",
            name,
            len(pka_data),
        )
        return None

    prot_states = []
    for i in range(0, len(pka_data), 3):
        try:
            site = pka_data[i]
            mean = float(pka_data[i + 1])
            std_base = float(pka_data[i + 2])
            logger.trace(
                "pKa data of site {}: Mean = {} ; Stdev = {}", site, mean, std_base
            )
            std = std_base * pka_std_range
            protonation_state = define_protonation_state(mean, std, min_ph, max_ph)
            logger.trace("Calculated protonation state of {}", protonation_state)
            prot_states.append([site, protonation_state])

        except (ValueError, IndexError) as e:
            logger.warning(
                "Error parsing pKa data for {} at position {}: {}", name, i, str(e)
            )
            return None

    return SubstructureData(
        name=name, smart=smart, mol=mol, prot_states_for_pH=prot_states
    )


def define_protonation_state(
    mean: float, std: float, min_ph: float, max_ph: float
) -> str:
    """Determine protonation state based on pKa and pH range.

    Args:
        mean: Mean pKa value
        std: Standard deviation (precision)
        min_ph: Minimum pH of range
        max_ph: Maximum pH of range

    Returns:
        Protonation state: "PROTONATED", "DEPROTONATED", or "BOTH"
    """
    min_pka = mean - std
    max_pka = mean + std

    # Check if pKa range overlaps with pH range
    if min_pka <= max_ph and min_ph <= max_pka:
        return "BOTH"
    elif mean > max_ph:
        return "PROTONATED"
    else:
        return "DEPROTONATED"


def get_prot_sites_and_target_states(
    smi: str, subs: list[SubstructureData]
) -> tuple[list[ProtonationSite], Mol | None]:
    """Find protonation sites and their target states for a molecule.

    Args:
        smi: SMILES string
        subs: List of substructure data

    Returns:
        Tuple of (protonation sites list, molecule object used for indexing)
    """
    logger.debug("Finding protonation sites for {}", smi)

    # Convert SMILES to mol object
    mol_used_to_idx_sites = smiles_to_mol(smi)
    if mol_used_to_idx_sites is None:
        return [], None

    # Add hydrogens
    try:
        logger.debug("Adding hydrogens")
        mol_with_hs = Chem.AddHs(mol_used_to_idx_sites)
        if mol_with_hs is None:
            logger.warning("Failed to add hydrogens to molecule: {}", smi)
            return [], None
        logger.trace("After adding hydrogens: {}", Chem.MolToSmiles(mol_with_hs))
        mol_used_to_idx_sites = mol_with_hs
    except Exception as e:
        logger.warning("Error adding hydrogens to molecule {} : {}", smi, str(e))
        return [], None

    # Unprotect molecule for substructure matching
    try:
        protect.unprotect_molecule(mol_used_to_idx_sites)
    except Exception as e:
        logger.warning("Error unprotecting molecule {} : {}", smi, str(e))
        return [], None

    protonation_sites = []
    matches_found = 0

    for sub_data in subs:
        if sub_data.mol is None:
            continue

        try:
            if mol_used_to_idx_sites.HasSubstructMatch(sub_data.mol):
                matches = protect.get_unprotected_matches(
                    mol_used_to_idx_sites, sub_data.mol
                )
                matches_found += len(matches)

                for match in matches:
                    for site_info in sub_data.prot_states_for_pH:
                        try:
                            site_idx = int(site_info[0])
                            target_state = site_info[1]

                            if site_idx >= len(match):
                                logger.warning(
                                    "Site index {} out of range for match in {}",
                                    site_idx,
                                    sub_data.name,
                                )
                                continue

                            prot_site = ProtonationSite(
                                atom_idx=match[site_idx],
                                target_state=target_state,
                                site_name=sub_data.name,
                            )

                            # Avoid duplicates
                            if not any(
                                ps.atom_idx == prot_site.atom_idx
                                and ps.target_state == prot_site.target_state
                                for ps in protonation_sites
                            ):
                                logger.debug(
                                    "Found {} at index {} with {} target",
                                    prot_site.site_name,
                                    prot_site.atom_idx,
                                    prot_site.target_state,
                                )
                                protonation_sites.append(prot_site)

                        except (ValueError, IndexError) as e:
                            logger.warning(
                                "Error processing site for {}: {}",
                                sub_data.name,
                                str(e),
                            )

                    # Protect this match to prevent overlapping matches
                    try:
                        protect.protect_molecule(mol_used_to_idx_sites, match)
                    except Exception as e:
                        logger.warning(
                            "Error protecting match for {}: {}", sub_data.name, str(e)
                        )

        except Exception as e:
            logger.warning(
                "Error matching substructure {} to '{}': {}", sub_data.name, smi, str(e)
            )

    logger.info("Found {} protonation site(s)", len(protonation_sites))

    return protonation_sites, mol_used_to_idx_sites


def protonate_site(mols: list[Mol], site: ProtonationSite) -> list[Mol]:
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

    logger.debug(
        "Protonating site {} ({}) with target state: {}",
        site.atom_idx,
        site.site_name,
        site.target_state,
    )

    # Define charge states for each protonation state
    state_to_charges = {"DEPROTONATED": [-1], "PROTONATED": [0], "BOTH": [-1, 0]}

    if site.target_state not in state_to_charges:
        logger.error("Invalid target protonation state: {}", site.target_state)
        return mols

    charges = state_to_charges[site.target_state]

    try:
        output_mols = set_protonation_charge(
            mols, site.atom_idx, charges, site.site_name
        )
        logger.debug("Generated {} protonated variants", len(output_mols))
        return output_mols
    except Exception as e:
        logger.error("Error protonating site {}: {}", site.atom_idx, str(e))
        return mols


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
    if not mols:
        return []

    output = []
    is_special_nitrogen = "*" in prot_site_name

    for charge in charges:
        nitrogen_charge = charge + 1

        # Special case for nitrogen moieties where acidic group is neutral
        if is_special_nitrogen:
            nitrogen_charge = nitrogen_charge - 1

        for mol in mols:
            try:
                processed_mol = _apply_charge_to_molecule(
                    mol, idx, charge, nitrogen_charge, prot_site_name
                )
                if processed_mol is not None:
                    output.append(processed_mol)
            except Exception as e:
                logger.warning(
                    "Error processing molecule with charge {}: {}", charge, str(e)
                )
                continue

    return output


def _apply_charge_to_molecule(
    mol: Mol, idx: int, charge: int, nitrogen_charge: int, prot_site_name: str
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


def validate_substructure_data(subs: list[SubstructureData]) -> dict[str, Any]:
    """Validate loaded substructure data and return statistics.

    Args:
        subs: List of substructure data

    Returns:
        Dictionary with validation statistics
    """
    stats: dict[str, int | dict[str, int]] = {
        "total_substructures": len(subs),
        "valid_mols": 0,
        "invalid_mols": 0,
        "total_sites": 0,
        "state_distribution": {"PROTONATED": 0, "DEPROTONATED": 0, "BOTH": 0},
    }

    for sub in subs:
        if sub.mol is not None:
            stats["valid_mols"] += 1
        else:
            stats["invalid_mols"] += 1

        stats["total_sites"] += len(sub.prot_states_for_pH)

        for site_info in sub.prot_states_for_pH:
            state = site_info[1]
            if state in stats["state_distribution"]:
                stats["state_distribution"][state] += 1

    return stats
