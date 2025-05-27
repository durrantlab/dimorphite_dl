"""
Detect protonation sites
"""

from typing import Iterable

import copy

from loguru import logger
from rdkit import Chem

from dimorphite_dl.mol import smiles_to_mol
from dimorphite_dl.protonate import PKaData
from dimorphite_dl.protonate.site import ProtonationSite


class ProtonationSiteDetector:
    @classmethod
    def find(cls, smiles: str) -> tuple[Chem.Mol | None, list[ProtonationSite]]:
        mol = cls._prep_mol(smiles)
        if mol is None:
            logger.warning("RDKit Mol conversion failed for {}", smiles)
            return None, []
        protonation_sites = cls._find_substruct_matches(mol)
        return mol, protonation_sites

    @classmethod
    def _get_mol(cls, smiles: str) -> Chem.Mol | None:
        logger.debug("Processing {}", smiles)
        mol = smiles_to_mol(smiles)
        mol = cls._add_hs(mol)

    @staticmethod
    def _add_hs(mol: Chem.Mol) -> Chem.Mol | None:
        logger.debug("Adding hydrogens")
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
    def _unprotect_atoms(mol):
        """Sets the protected property on all atoms to 0."""
        logger.trace("Unprotecting each atom")
        for atom in mol.GetAtoms():
            atom.SetProp("_protected", "0")

    @classmethod
    def _prep_mol(cls, smiles: str) -> Chem.Mol | None:
        mol = cls._get_mol(smiles)
        mol = cls._add_hs(mol)
        mol = cls._unprotect_atoms(mol)
        return mol

    @staticmethod
    def _protect_match(mol: Chem.Mol, matches: list[int]):
        """Given a 'match', a list of molecules idx's, we set the protected status
        of each atom to 1."""
        logger.trace("Protecting atom(s)")
        for idx in matches:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp("_protected", "1")

    @staticmethod
    def is_match_unprotected(mol: Chem.Mol, match: Iterable[int]) -> bool:
        """Checks a molecule to see if the substructure match contains any
        protected atoms.

        Args:
            mol: RDKit Mol to check.
            match: Match to check.

        Returns:
            Whether the match is present or not.
        """
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetProp("_protected")
            if protected == "1":
                return False
        return True

    @classmethod
    def get_unprotected_matches(cls, mol: Chem.Mol, substruct: Chem.Mol) -> list[int]:
        """Finds substructure matches with atoms that have not been protected.

        Args:
            mol: Prepped RDKit Mol object from SMILES.
            substruct: Substructure to search for in `mol`

        Returns:
            Atom indices that match the `substruct` query.
        """
        matches = mol.GetSubstructMatches(substruct)
        logger.debug("Found {} substructure match(es)", len(matches))
        unprotected_matches = []
        for match in matches:
            if cls.is_match_unprotected(mol, match):
                unprotected_matches.append(match)
        logger.debug(
            "{}/{} matches were unprotected", len(unprotected_matches), len(matches)
        )
        return unprotected_matches

    @classmethod
    def _get_substruct_matches(cls, mol: Chem.Mol, sub_mol: Chem.Mol) -> list[int]:
        if mol.HasSubstructMatch(sub_mol):
            matches = cls.get_unprotected_matches(mol, sub_mol)
        else:
            matches = []
        return matches

    @classmethod
    def _find_substruct_matches(cls, mol: Chem.Mol) -> list[ProtonationSite]:
        site_matches = []
        matches_found = 0

        for sub_data in PKaData.get_substructures():
            if sub_data.mol is None:
                continue

            matches = cls._get_substruct_matches(mol, sub_data.mol)
            if len(matches) == 0:
                continue
            else:
                matches_found += len(matches)

            for match in matches:
                # Maybe need to check for duplicates?
                site_matches.append(
                    ProtonationSite(
                        idx_atom=match, substructure=copy.deepcopy(sub_data)
                    )
                )
            _mol = cls._protect_match(mol, matches)
            if _mol is not None:
                mol = _mol
            else:
                logger.warning("Could not protect match at index {}", matches)
        logger.info("Found {} protonation site(s)", matches_found)
        return site_matches
