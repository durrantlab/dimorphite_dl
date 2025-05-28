"""
This module provides functionality to detect protonation sites in molecules
using substructure matching with comprehensive error handling and validation.
"""

from typing import Generator

import copy

from loguru import logger
from rdkit import Chem

from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.protonate.data import PKaData
from dimorphite_dl.protonate.site import ProtonationSite, SubstructureDatum


class ProtonationSiteDetectionError(Exception):
    """Raised when protonation site detection encounters an error."""

    pass


class ProtonationSiteDetector:
    """
    Robust detector for finding protonation sites in molecules.

    Uses substructure matching to identify potential protonation sites
    based on known pKa patterns. Handles atom protection to prevent
    overlapping matches.
    """

    def __init__(self, validate_sites: bool = True, max_sites_per_molecule: int = 50):
        """
        Initialize the detector with explicit configuration.

        Args:
            validate_sites: Whether to validate detected sites (explicit, not default)
            max_sites_per_molecule: Maximum sites to detect per molecule (bounded)
        """
        assert isinstance(validate_sites, bool)
        assert isinstance(max_sites_per_molecule, int)
        assert max_sites_per_molecule > 0
        assert max_sites_per_molecule <= 1000  # Reasonable upper bound

        self.validate_sites = validate_sites
        self.max_sites_per_molecule = max_sites_per_molecule
        self.pka_data = PKaData()

        # Initialize statistics - all counters start at zero
        self._stats_molecules_processed = 0
        self._stats_sites_found = 0
        self._stats_sites_validated = 0
        self._stats_sites_rejected = 0
        self._stats_substructures_matched = 0

    def find_sites(
        self, mol_record: MoleculeRecord
    ) -> tuple[MoleculeRecord, list[ProtonationSite]]:
        """
        Find protonation sites in a molecule. This is the main entry point.

        Args:
            mol_record: MoleculeRecord to analyze

        Returns:
            Tuple of (updated_mol_record, list_of_protonation_sites)

        Raises:
            ProtonationSiteDetectionError: If detection fails critically
        """
        assert isinstance(mol_record, MoleculeRecord)
        assert mol_record.smiles  # SMILES cannot be empty

        logger.debug("Finding protonation sites for '{}'", mol_record.smiles)
        self._stats_molecules_processed += 1

        try:
            prepared_mol = mol_record.prepare_for_protonation()
            if prepared_mol is None:
                logger.warning("Failed to prepare molecule: '{}'", mol_record.smiles)
                return mol_record, []

            mol_record.mol = prepared_mol
            protonation_sites = self._detect_all_sites_in_molecule(prepared_mol)

            if self.validate_sites:
                protonation_sites = self._validate_detected_sites(
                    protonation_sites, mol_record
                )

            sites_count = len(protonation_sites)
            self._stats_sites_found += sites_count

            logger.info(
                "Found {} protonation site(s) for '{}'", sites_count, mol_record.smiles
            )
            return mol_record, protonation_sites

        except Exception as error:
            logger.error(
                "Critical error detecting sites for '{}': {}",
                mol_record.smiles,
                str(error),
            )
            raise ProtonationSiteDetectionError(f"Detection failed: {error}") from error

    def _detect_all_sites_in_molecule(self, mol: Chem.Mol) -> list[ProtonationSite]:
        """
        Detect all protonation sites in the prepared molecule.

        Args:
            mol: Prepared RDKit mol object

        Returns:
            List of detected protonation sites
        """
        assert mol is not None
        assert mol.GetNumAtoms() > 0

        protonation_sites = []
        total_matches_found = 0

        for substructure_data in self._iterate_available_substructures():
            if substructure_data.mol is None:
                logger.debug(
                    "Skipping substructure '{}' - no mol object", substructure_data.name
                )
                continue

            matches = self._find_unprotected_matches_for_substructure(
                mol, substructure_data
            )

            if len(matches) == 0:
                continue

            match_count = len(matches)
            total_matches_found += match_count
            self._stats_substructures_matched += 1

            logger.debug(
                "Found {} match(es) for substructure '{}'",
                match_count,
                substructure_data.name,
            )

            sites_created = self._create_sites_from_matches(
                matches, substructure_data, protonation_sites
            )
            if sites_created == 0:
                break  # Hit the limit, stop processing

            mol = self._protect_matched_atoms_in_molecule(mol, matches)

            if len(protonation_sites) >= self.max_sites_per_molecule:
                break

        logger.debug(
            "Total matches found: {}, sites created: {}",
            total_matches_found,
            len(protonation_sites),
        )
        assert len(protonation_sites) <= self.max_sites_per_molecule

        return protonation_sites

    def _iterate_available_substructures(
        self,
    ) -> Generator[SubstructureDatum, None, None]:
        """
        Get available substructure patterns for matching.

        Yields:
            SubstructureDatum objects for pattern matching
        """
        try:
            substructure_count = 0
            for substructure_data in self.pka_data.get_substructures():
                assert isinstance(substructure_data, SubstructureDatum)
                substructure_count += 1
                yield substructure_data

            logger.trace("Iterated over {} substructures", substructure_count)

        except Exception as error:
            logger.error("Error loading substructure data: {}", str(error))
            raise ProtonationSiteDetectionError(
                f"Failed to load substructure data: {error}"
            ) from error

    def _find_unprotected_matches_for_substructure(
        self, mol: Chem.Mol, substructure_data: SubstructureDatum
    ) -> list[tuple[int, ...]]:
        """
        Find unprotected matches for a specific substructure pattern.

        Args:
            mol: RDKit mol object to search in
            substructure_data: Substructure pattern to match

        Returns:
            List of tuples containing atom indices for unprotected matches
        """
        assert mol is not None
        assert substructure_data is not None
        assert substructure_data.mol is not None

        try:
            has_substructure = mol.HasSubstructMatch(substructure_data.mol)
            if not has_substructure:
                return []

            all_matches = list(mol.GetSubstructMatches(substructure_data.mol))
            total_matches = len(all_matches)
            logger.debug(
                "Found {} total match(es) for '{}'",
                total_matches,
                substructure_data.name,
            )

            unprotected_matches = self._filter_matches_by_protection_status(
                mol, all_matches
            )
            unprotected_count = len(unprotected_matches)

            logger.debug(
                "{}/{} matches were unprotected for '{}'",
                unprotected_count,
                total_matches,
                substructure_data.name,
            )

            return unprotected_matches

        except Exception as error:
            logger.warning(
                "Error finding matches for substructure '{}': {}",
                substructure_data.name,
                str(error),
            )
            return []

    def _filter_matches_by_protection_status(
        self, mol: Chem.Mol, all_matches: list[tuple[int, ...]]
    ) -> list[tuple[int, ...]]:
        """
        Filter matches to only include those with unprotected atoms.

        Args:
            mol: RDKit mol object
            all_matches: List of all matches to filter

        Returns:
            List of matches where all atoms are unprotected
        """
        assert mol is not None
        assert isinstance(all_matches, list)

        unprotected_matches = []
        atom_count = mol.GetNumAtoms()

        for match in all_matches:
            assert isinstance(match, tuple)

            # Validate atom indices are within bounds
            for atom_index in match:
                assert isinstance(atom_index, int)
                assert 0 <= atom_index < atom_count

            if self._are_all_atoms_in_match_unprotected(mol, match):
                unprotected_matches.append(match)

        return unprotected_matches

    def _are_all_atoms_in_match_unprotected(
        self, mol: Chem.Mol, match: tuple[int, ...]
    ) -> bool:
        """
        Check if all atoms in a match are unprotected.

        Args:
            mol: RDKit mol object
            match: Tuple of atom indices to check

        Returns:
            True if all atoms in match are unprotected, False otherwise
        """
        assert mol is not None
        assert isinstance(match, tuple)
        assert len(match) > 0  # Match cannot be empty

        try:
            for atom_index in match:
                assert isinstance(atom_index, int)
                if MoleculeRecord.is_atom_protected(mol, atom_index):
                    return False
            return True

        except Exception as error:
            logger.debug(
                "Error checking protection for match {}: {}", match, str(error)
            )
            return False

    def _create_sites_from_matches(
        self,
        matches: list[tuple[int, ...]],
        substructure_data: SubstructureDatum,
        existing_sites: list[ProtonationSite],
    ) -> int:
        """
        Create ProtonationSite objects from matches.

        Args:
            matches: List of atom index tuples
            substructure_data: Substructure information
            existing_sites: Current list of sites (modified in place)

        Returns:
            Number of sites actually created (may be limited by max_sites)
        """
        assert isinstance(matches, list)
        assert isinstance(substructure_data, SubstructureDatum)
        assert isinstance(existing_sites, list)

        sites_created = 0

        for match_indices in matches:
            current_site_count = len(existing_sites)
            if current_site_count >= self.max_sites_per_molecule:
                logger.warning(
                    "Reached maximum sites limit ({}) for molecule",
                    self.max_sites_per_molecule,
                )
                break

            site = ProtonationSite(
                idx_atom=match_indices, substructure=copy.deepcopy(substructure_data)
            )
            existing_sites.append(site)
            sites_created += 1

        assert sites_created <= len(matches)
        assert len(existing_sites) <= self.max_sites_per_molecule

        return sites_created

    def _protect_matched_atoms_in_molecule(
        self, mol: Chem.Mol, matches: list[tuple[int, ...]]
    ) -> Chem.Mol:
        """
        Protect all atoms involved in matches to prevent overlap.

        Args:
            mol: RDKit mol object
            matches: List of matches whose atoms should be protected

        Returns:
            Same mol object with matched atoms protected
        """
        assert mol is not None
        assert isinstance(matches, list)

        for match in matches:
            assert isinstance(match, tuple)
            atom_indices = list(match)
            mol = MoleculeRecord.protect_atoms(mol, atom_indices)

        return mol

    def _validate_detected_sites(
        self, sites: list[ProtonationSite], mol_record: MoleculeRecord
    ) -> list[ProtonationSite]:
        """
        Validate all detected protonation sites.

        Args:
            sites: List of protonation sites to validate
            mol_record: Original molecule record for context

        Returns:
            List of validated protonation sites
        """
        assert isinstance(sites, list)
        assert isinstance(mol_record, MoleculeRecord)

        if len(sites) == 0:
            return sites

        logger.debug("Validating {} protonation sites", len(sites))
        validated_sites = []

        for site_index, site in enumerate(sites):
            assert isinstance(site, ProtonationSite)

            try:
                if self._is_single_site_valid(site, mol_record):
                    validated_sites.append(site)
                    self._stats_sites_validated += 1
                else:
                    self._stats_sites_rejected += 1
                    logger.debug(
                        "Rejected site {} due to validation failure", site_index
                    )

            except Exception as error:
                logger.warning("Error validating site {}: {}", site_index, str(error))
                self._stats_sites_rejected += 1

        validated_count = len(validated_sites)
        original_count = len(sites)
        logger.debug("Validated {}/{} sites", validated_count, original_count)

        assert validated_count <= original_count
        return validated_sites

    def _is_single_site_valid(
        self, site: ProtonationSite, mol_record: MoleculeRecord
    ) -> bool:
        """
        Validate a single protonation site.

        Args:
            site: ProtonationSite to validate
            mol_record: Original molecule record for bounds checking

        Returns:
            True if site is valid, False otherwise
        """
        assert isinstance(site, ProtonationSite)
        assert isinstance(mol_record, MoleculeRecord)

        mol = mol_record.mol
        if mol is None:
            logger.debug("Site validation failed: no mol object")
            return False

        if not self._are_atom_indices_valid_for_molecule(site.idx_atom, mol):
            return False

        if not self._is_substructure_data_valid(site.substructure):
            return False

        return True

    def _are_atom_indices_valid_for_molecule(self, idx_atom, mol: Chem.Mol) -> bool:
        """
        Check if atom indices are valid for the molecule.

        Args:
            idx_atom: Single index or tuple of indices
            mol: RDKit mol object for bounds checking

        Returns:
            True if all indices are valid, False otherwise
        """
        assert mol is not None

        atom_count = mol.GetNumAtoms()
        assert atom_count > 0

        # Handle both single atom index and tuple of indices
        if isinstance(idx_atom, (list, tuple)):
            indices_to_check = idx_atom
        else:
            indices_to_check = [idx_atom]

        for atom_index in indices_to_check:
            if not isinstance(atom_index, int):
                logger.debug("Invalid atom index type: {}", type(atom_index))
                return False
            if atom_index < 0:
                logger.debug("Negative atom index: {}", atom_index)
                return False
            if atom_index >= atom_count:
                logger.debug(
                    "Atom index {} out of range (molecule has {} atoms)",
                    atom_index,
                    atom_count,
                )
                return False

        return True

    def _is_substructure_data_valid(self, substructure: SubstructureDatum) -> bool:
        """
        Check if substructure data is valid.

        Args:
            substructure: SubstructureDatum to validate

        Returns:
            True if substructure data is valid, False otherwise
        """
        if substructure is None:
            logger.debug("Site has null substructure")
            return False

        if not substructure.pkas:
            logger.debug("Site has no pKa data")
            return False

        if len(substructure.pkas) == 0:
            logger.debug("Site has empty pKa list")
            return False

        return True

    def get_stats(self) -> dict[str, int]:
        """
        Get detection statistics.

        Returns:
            Dictionary of detection statistics
        """
        return {
            "molecules_processed": self._stats_molecules_processed,
            "sites_found": self._stats_sites_found,
            "sites_validated": self._stats_sites_validated,
            "sites_rejected": self._stats_sites_rejected,
            "substructures_matched": self._stats_substructures_matched,
        }

    def reset_stats(self) -> None:
        """
        Reset all detection statistics to zero.

        """
        self._stats_molecules_processed = 0
        self._stats_sites_found = 0
        self._stats_sites_validated = 0
        self._stats_sites_rejected = 0
        self._stats_substructures_matched = 0


def canonicalize_smiles_list(
    mols: list[Chem.Mol], original_smiles: str = ""
) -> list[str]:
    """
    Generate canonical SMILES from molecule objects.

    Args:
        mols: List of RDKit mol objects to convert
        original_smiles: Original SMILES for logging context

    Returns:
        List of unique canonical SMILES strings
    """
    assert isinstance(mols, list)
    assert isinstance(original_smiles, str)

    if len(mols) == 0:
        return []

    logger.debug("Generating canonical SMILES for {} molecules", len(mols))

    try:
        unique_smiles = set()
        valid_mol_count = 0

        for mol in mols:
            if mol is None:
                continue

            valid_mol_count += 1
            canonical_smiles = _generate_single_canonical_smiles(mol)
            if canonical_smiles is not None:
                unique_smiles.add(canonical_smiles)

        smiles_list = list(unique_smiles)
        unique_count = len(smiles_list)

        context_msg = f" for '{original_smiles}'" if original_smiles else ""
        logger.debug(
            "Generated {} unique canonical SMILES from {} valid molecules{}",
            unique_count,
            valid_mol_count,
            context_msg,
        )

        return smiles_list

    except Exception as error:
        context_msg = f" for '{original_smiles}'" if original_smiles else ""
        logger.warning(
            "Error in canonical SMILES generation{}: {}", context_msg, str(error)
        )
        return []


def _generate_single_canonical_smiles(mol: Chem.Mol) -> str | None:
    """
    Generate canonical SMILES for a single molecule.

    Args:
        mol: RDKit mol object

    Returns:
        Canonical SMILES string or None if generation failed
    """
    assert mol is not None

    try:
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        if canonical and len(canonical) > 0:
            return canonical
        else:
            logger.debug("Generated empty SMILES string")
            return None

    except Exception as error:
        logger.debug("Error generating canonical SMILES: {}", str(error))
        return None


# Convenience functions for backward compatibility
def find(mol_record: MoleculeRecord) -> tuple[MoleculeRecord, list[ProtonationSite]]:
    """
    Convenience function for finding protonation sites with default settings.

    Args:
        mol_record: MoleculeRecord to analyze

    Returns:
        Tuple of (updated_mol_record, list_of_protonation_sites)
    """
    assert isinstance(mol_record, MoleculeRecord)

    detector = ProtonationSiteDetector(validate_sites=True, max_sites_per_molecule=50)
    return detector.find_sites(mol_record)


def find_with_options(
    mol_record: MoleculeRecord, validate_sites: bool = True, max_sites: int = 50
) -> tuple[MoleculeRecord, list[ProtonationSite]]:
    """
    Convenience function for finding protonation sites with explicit options.

    Args:
        mol_record: MoleculeRecord to analyze
        validate_sites: Whether to validate detected sites (explicit)
        max_sites: Maximum sites to detect per molecule (bounded)

    Returns:
        Tuple of (updated_mol_record, list_of_protonation_sites)
    """
    assert isinstance(mol_record, MoleculeRecord)
    assert isinstance(validate_sites, bool)
    assert isinstance(max_sites, int)
    assert max_sites > 0

    detector = ProtonationSiteDetector(
        validate_sites=validate_sites, max_sites_per_molecule=max_sites
    )
    return detector.find_sites(mol_record)
