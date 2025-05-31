"""
Enhanced MoleculeRecord class for handling SMILES strings and RDKit mol objects.

This class consolidates all molecule-related operations including:
- SMILES validation and conversion
- RDKit mol object management
- Neutralization
- Hydrogen handling
- Error handling and logging
"""

from typing import Any

import copy
import os
import sys

from loguru import logger
from rdkit import Chem

from dimorphite_dl.neutralize import MoleculeNeutralizer


class MoleculeRecord:
    """
    Enhanced class for managing SMILES strings and RDKit mol objects.

    Handles all molecule-related operations including validation, conversion,
    neutralization, and hydrogen management.
    """

    def __init__(self, smiles: str, identifier: str = "") -> None:
        """
        Initialize a MoleculeRecord.

        Args:
            smiles: SMILES string representation of the molecule
            identifier: Optional unique identifier for the molecule

        Raises:
            ValueError: If smiles is not a valid string
        """
        if not isinstance(smiles, str):
            raise ValueError(f"SMILES must be a string, got {type(smiles)}")
        if not isinstance(identifier, str):
            raise ValueError(f"Identifier must be a string, got {type(identifier)}")

        smiles = smiles.strip()
        if not smiles:
            raise ValueError("SMILES string cannot be empty")

        self.smiles_original = smiles
        """Original SMILES used to initialize this MoleculeRecord"""

        self.identifier = identifier
        """Unique identifier for molecule"""

        self.smiles = smiles
        """Current SMILES after any processing"""

        self._mol: Chem.Mol | None = None
        """Cached RDKit mol object"""

        self._mol_with_hs: Chem.Mol | None = None
        """Cached RDKit mol object with explicit hydrogens"""

        self._neutralizer: MoleculeNeutralizer | None = None
        """Cached neutralizer instance"""

    @property
    def mol(self) -> Chem.Mol | None:
        """Get the RDKit mol object, creating it if necessary."""
        if self._mol is None:
            self._mol = self.to_mol()
        return self._mol

    @mol.setter
    def mol(self, value: Chem.Mol | None) -> None:
        """Set the RDKit mol object and clear dependent caches."""
        self._mol = value
        self._mol_with_hs = None  # Clear dependent cache

    def to_mol(self) -> Chem.Mol | None:
        """
        Convert current SMILES to a RDKit Mol object.

        Returns:
            RDKit Mol object or None if conversion fails
        """
        conversion_info = self.to_mol_silenced(self.smiles)

        if conversion_info["mol"] is None:
            error_msg = conversion_info["stderr_content"].strip()
            if error_msg:
                logger.warning(
                    "RDKit failed to parse SMILES '{}'. RDKit error: {}",
                    self.smiles,
                    error_msg,
                )
            else:
                logger.warning(
                    "RDKit failed to parse SMILES '{}' (no specific error message)",
                    self.smiles,
                )
            return None

        mol = conversion_info["mol"]
        logger.trace("SMILES after conversion: {}", Chem.MolToSmiles(mol))
        return mol

    def to_mol_with_hs(self) -> Chem.Mol | None:
        """
        Get RDKit mol object with explicit hydrogens.

        Returns:
            RDKit Mol object with explicit hydrogens or None if conversion fails
        """
        if self._mol_with_hs is None:
            base_mol = self.mol
            if base_mol is not None:
                self._mol_with_hs = self.add_hydrogens(base_mol)
        return self._mol_with_hs

    def refresh_mol_from_smiles(self) -> bool:
        """
        Refresh the mol object from current SMILES string.

        Returns:
            True if successful, False otherwise
        """
        self._mol = None
        self._mol_with_hs = None
        new_mol = self.to_mol()
        return new_mol is not None

    def update_smiles_from_mol(self, mol: Chem.Mol | None = None) -> bool:
        """
        Update SMILES string from RDKit mol object.

        Args:
            mol: Optional mol object to use. If None, uses self.mol

        Returns:
            True if successful, False otherwise
        """
        if mol is None:
            mol = self.mol

        if mol is None:
            logger.warning("Cannot update SMILES: no valid mol object available")
            return False

        try:
            new_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if new_smiles:
                self.smiles = new_smiles
                # Clear mol cache since we're updating from external mol
                if mol is not self._mol:
                    self._mol = mol
                    self._mol_with_hs = None
                return True
        except Exception as e:
            logger.warning("Error generating SMILES from mol: {}", str(e))

        return False

    def get_canonical_smiles(self, isomeric: bool = True) -> str | None:
        """
        Get canonical SMILES representation.

        Args:
            isomeric: Whether to include stereochemistry information

        Returns:
            Canonical SMILES string or None if conversion fails
        """
        mol = self.mol
        if mol is None:
            return None

        try:
            return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=True)
        except Exception as e:
            logger.warning("Error generating canonical SMILES: {}", str(e))
            return None

    def make_canonical(self, isomeric: bool = True) -> None:
        smiles = self.get_canonical_smiles(isomeric=isomeric)
        if smiles is not None:
            self._update_smiles(smiles)

    def is_valid(self) -> bool:
        """
        Check if the current SMILES represents a valid molecule.

        Returns:
            True if valid, False otherwise
        """
        return self.mol is not None

    def _update_smiles(self, smiles: str) -> None:
        self.smiles = smiles
        # Clear caches since SMILES changed
        self._mol = None
        self._mol_with_hs = None

    def get_neutralized(self, smiles: str) -> str:
        if self._neutralizer is None:
            self._neutralizer = MoleculeNeutralizer()

        neutralized_smiles = self._neutralizer.neutralize_smiles(smiles)
        if neutralized_smiles is not None:
            logger.debug("Successfully neutralized molecule")
            return neutralized_smiles
        raise RuntimeError("Issue neutralizing SMILES")

    def neutralize(self):
        """
        Neutralize the molecule using the neutralizer.

        Returns:
            True if neutralization was successful, False otherwise
        """
        smiles_neutralized = self.get_neutralized(self.smiles)
        self._update_smiles(smiles_neutralized)

    @staticmethod
    def add_hydrogens(mol: Chem.Mol) -> Chem.Mol | None:
        """
        Add explicit hydrogens to a molecule.

        Args:
            mol: RDKit mol object

        Returns:
            Mol object with explicit hydrogens or None if failed
        """
        if mol is None:
            return None

        logger.debug("Adding hydrogens to molecule")
        try:
            mol_with_hs = Chem.AddHs(mol)
            if mol_with_hs is None:
                logger.warning("Failed to add hydrogens to molecule")
                return None
            logger.trace("After adding hydrogens: {}", Chem.MolToSmiles(mol_with_hs))
            return mol_with_hs
        except Exception as e:
            logger.warning("Error adding hydrogens to molecule: {}", str(e))
            return None

    @staticmethod
    def remove_hydrogens(mol: Chem.Mol) -> Chem.Mol | None:
        """
        Remove explicit hydrogens from a molecule.

        Args:
            mol: RDKit mol object

        Returns:
            Mol object without explicit hydrogens or None if failed
        """
        if mol is None:
            return None

        logger.debug("Removing hydrogens from molecule")
        try:
            mol_no_hs = Chem.RemoveHs(mol)
            if mol_no_hs is None:
                logger.warning("Failed to remove hydrogens from molecule")
                return None
            return mol_no_hs
        except Exception as e:
            logger.warning("Error removing hydrogens from molecule: {}", str(e))
            return None

    @staticmethod
    def unprotect_atoms(mol: Chem.Mol) -> Chem.Mol:
        """
        Set the protected property on all atoms to 0.

        Args:
            mol: RDKit mol object to unprotect

        Returns:
            The same mol object (modified in place)
        """
        logger.trace("Unprotecting each atom")
        for atom in mol.GetAtoms():
            atom.SetProp("_protected", "0")
        return mol

    @staticmethod
    def protect_atoms(mol: Chem.Mol, atom_indices: list[int]) -> Chem.Mol:
        """
        Set the protected property on specified atoms to 1.

        Args:
            mol: RDKit mol object
            atom_indices: List of atom indices to protect

        Returns:
            The same mol object (modified in place)
        """
        logger.trace("Protecting atom(s): {}", atom_indices)
        for idx in atom_indices:
            try:
                atom = mol.GetAtomWithIdx(idx)
                atom.SetProp("_protected", "1")
            except Exception as e:
                logger.warning("Could not protect atom at index {}: {}", idx, str(e))
        return mol

    @staticmethod
    def is_atom_protected(mol: Chem.Mol, atom_idx: int) -> bool:
        """
        Check if an atom is protected.

        Args:
            mol: RDKit mol object
            atom_idx: Atom index to check

        Returns:
            True if atom is protected, False otherwise
        """
        try:
            atom = mol.GetAtomWithIdx(atom_idx)
            protected = atom.GetProp("_protected")
            return protected == "1"
        except Exception:
            return False

    def process_azides(self) -> None:
        """
        Process azide patterns in SMILES string.

        Args:
            smiles: Input SMILES string

        Returns:
            SMILES string with processed azides
        """
        smiles_working = self.smiles
        if "N=N=N" in smiles_working or "NN#N" in smiles_working:
            logger.info("Attempting to fix azide patterns in: '{}'", smiles_working)
            smiles_working = smiles_working.replace("N=N=N", "N=[N+]=N")
            smiles_working = smiles_working.replace("NN#N", "N=[N+]=N")
            if smiles_working != self.smiles:
                logger.info(
                    "Modified SMILES: '{}' -> '{}'", self.smiles, smiles_working
                )
        self._update_smiles(smiles_working)

    def prepare_for_protonation(self) -> Chem.Mol:
        """
        Prepare molecule for protonation site detection.

        Returns:
            Prepared RDKit mol object or None if preparation fails
        """
        logger.info("Preparing molecule for analysis")

        self.process_azides()
        self.neutralize()

        base_mol = self.to_mol()
        if base_mol is None:
            raise RuntimeError("Could not convert SMILES to RDKit Mol")

        mol_with_hydrogens = self.add_hydrogens(base_mol)
        if mol_with_hydrogens is None:
            raise RuntimeError("Could not add Hydrogens to Mol")

        prepared_mol = self.unprotect_atoms(mol_with_hydrogens)

        atom_count = prepared_mol.GetNumAtoms()
        logger.trace("Molecule prepared with {} atoms", atom_count)
        assert atom_count > 0  # Molecule must have at least one atom

        self._update_smiles(Chem.MolToSmiles(prepared_mol))
        return prepared_mol

    @staticmethod
    def to_mol_silenced(smiles: str) -> dict[str, Any]:
        """
        Capture RDKit stderr output and return mol object with error messages.

        Args:
            smiles: SMILES string to convert

        Returns:
            Dictionary with 'mol' (RDKit Mol or None) and 'stderr_content' (string)
        """
        logger.debug("Converting SMILES to RDKit mol: {}", smiles)

        # Set up stderr capture
        stderr_fileno = sys.stderr.fileno()
        stderr_save = os.dup(stderr_fileno)
        stderr_pipe = os.pipe()

        try:
            # Redirect stderr to pipe
            os.dup2(stderr_pipe[1], stderr_fileno)
            os.close(stderr_pipe[1])

            # Convert SMILES to mol (this may write to stderr)
            mol = Chem.MolFromSmiles(smiles)

            # Read captured stderr
            os.close(stderr_fileno)
            stderr_content = os.read(stderr_pipe[0], 1024).decode(
                "utf-8", errors="ignore"
            )

        except Exception as e:
            logger.error("Error during SMILES conversion: {}", str(e))
            mol = None
            stderr_content = f"Exception during conversion: {str(e)}"

        finally:
            # Restore stderr
            try:
                os.close(stderr_pipe[0])
            except Exception:
                pass
            try:
                os.dup2(stderr_save, stderr_fileno)
                os.close(stderr_save)
            except Exception:
                pass

        return {"mol": mol, "stderr_content": stderr_content}

    def copy(self) -> "MoleculeRecord":
        """
        Create a deep copy of this MoleculeRecord.

        Returns:
            New MoleculeRecord instance
        """
        new_record = MoleculeRecord(self.smiles, self.identifier)
        new_record.smiles_original = self.smiles_original

        # Deep copy mol objects if they exist
        if self._mol is not None:
            new_record._mol = copy.deepcopy(self._mol)
        if self._mol_with_hs is not None:
            new_record._mol_with_hs = copy.deepcopy(self._mol_with_hs)

        return new_record

    def get_atom_count(self) -> int:
        """
        Get the number of atoms in the molecule.

        Returns:
            Number of atoms, or 0 if mol is invalid
        """
        mol = self.mol
        return mol.GetNumAtoms() if mol is not None else 0

    def get_heavy_atom_count(self) -> int:
        """
        Get the number of heavy (non-hydrogen) atoms in the molecule.

        Returns:
            Number of heavy atoms, or 0 if mol is invalid
        """
        mol = self.mol
        return mol.GetNumHeavyAtoms() if mol is not None else 0

    def has_substructure(self, pattern: str | Chem.Mol) -> bool:
        """
        Check if molecule contains a specific substructure.

        Args:
            pattern: SMARTS string or RDKit mol object to search for

        Returns:
            True if substructure is found, False otherwise
        """
        mol = self.mol
        if mol is None:
            return False

        try:
            if isinstance(pattern, str):
                pattern_mol = Chem.MolFromSmarts(pattern)
                if pattern_mol is None:
                    logger.warning("Invalid SMARTS pattern: {}", pattern)
                    return False
            else:
                pattern_mol = pattern

            return mol.HasSubstructMatch(pattern_mol)
        except Exception as e:
            logger.warning("Error checking substructure: {}", str(e))
            return False

    def get_substructure_matches(
        self, pattern: str | Chem.Mol
    ) -> list[tuple[int, ...]]:
        """
        Get all matches of a substructure pattern.

        Args:
            pattern: SMARTS string or RDKit mol object to search for

        Returns:
            List of tuples containing atom indices for each match
        """
        mol = self.mol
        if mol is None:
            return []

        try:
            if isinstance(pattern, str):
                pattern_mol = Chem.MolFromSmarts(pattern)
                if pattern_mol is None:
                    logger.warning("Invalid SMARTS pattern: {}", pattern)
                    return []
            else:
                pattern_mol = pattern

            return list(mol.GetSubstructMatches(pattern_mol))
        except Exception as e:
            logger.warning("Error finding substructure matches: {}", str(e))
            return []

    def __str__(self) -> str:
        """String representation of the molecule."""
        if self.identifier:
            return f"MoleculeRecord('{self.smiles}', '{self.identifier}')"
        return f"MoleculeRecord('{self.smiles}')"

    def __repr__(self) -> str:
        """Detailed string representation of the molecule."""
        return self.__str__()

    def __eq__(self, other: object) -> bool:
        """Check equality based on canonical SMILES."""
        if not isinstance(other, MoleculeRecord):
            return False

        self_canonical = self.get_canonical_smiles()
        other_canonical = other.get_canonical_smiles()

        return (
            self_canonical is not None
            and other_canonical is not None
            and self_canonical == other_canonical
        )

    def __hash__(self) -> int:
        """Hash based on canonical SMILES."""
        canonical = self.get_canonical_smiles()
        return hash(canonical) if canonical is not None else hash(self.smiles)
