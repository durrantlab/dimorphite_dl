import os
import sys

from loguru import logger
from rdkit import Chem


class MoleculeRecord:
    def __init__(self, smiles: str, identifier: str = "") -> None:
        assert isinstance(smiles, str)
        assert isinstance(identifier, str)
        smiles = smiles.strip()

        self.smiles_original = smiles
        """Original SMILES used to initialize this `MoleculeRecord`"""

        self.identifier = identifier
        """Unique identifier for molecule"""

        self.smiles = smiles
        """Current SMILES after any processing."""

    def to_mol(self) -> Chem.Mol | None:
        """Converts current `self.smiles` to a RDKit Mol. Does no processing"""
        conversion_info = MoleculeRecord.to_mol_silenced(self.smiles)
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

        logger.trace(
            "SMILES after conversion: {}", Chem.MolToSmiles(conversion_info["mol"])
        )
        return conversion_info["mol"]

    @staticmethod
    def _process_azides(smiles: str) -> str:
        """Process azides in SMILES string"""
        smiles_working = smiles
        if "N=N=N" in smiles_working or "NN#N" in smiles_working:
            logger.info("Attempting to fix azide patterns in: '{}'", smiles_working)
            smiles_working = smiles_working.replace("N=N=N", "N=[N+]=N")
            smiles_working = smiles_working.replace("NN#N", "N=[N+]=N")
            if smiles_working != smiles:
                logger.info("Modified SMILES: '{}' -> '{}'", smiles, smiles_working)
        return smiles_working

    @staticmethod
    def to_mol_silenced(smiles: str) -> dict[str, Chem.Mol | str]:
        """Capture RDKit stderr output and return mol object with error messages.

        Args:
            smiles_str: SMILES string to convert

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

    @staticmethod
    def neutralize_smiles(smiles: str) -> str:
        """Neutralize SMILES"""
        logger.debug("Neutralizing {}", smiles)
        # Get the reaction data
        rxn_data = [
            [
                "[Ov1-1:1]",
                "[Ov2+0:1]-[H]",
            ],  # To handle O- bonded to only one atom (add hydrogen).
            [
                "[#7v4+1:1]-[H]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to a hydrogen (remove hydrogen).
            [
                "[Ov2-:1]",
                "[Ov2+0:1]",
            ],  # To handle O- bonded to two atoms. Should not be Negative.
            [
                "[#7v3+1:1]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to three atoms. Should not be positive.
            [
                "[#7v2-1:1]",
                "[#7+0:1]-[H]",
            ],  # To handle N- Bonded to two atoms. Add hydrogen.
            # ['[N:1]=[N+0:2]=[N:3]-[H]', '[N:1]=[N+1:2]=[N+0:3]-[H]'],  # To handle bad azide. Must be
            # protonated. (Now handled
            # elsewhere, before SMILES
            # converted to Mol object.)
            [
                "[H]-[N:1]-[N:2]#[N:3]",
                "[N:1]=[N+1:2]=[N:3]-[H]",
            ],  # To handle bad azide. R-N-N#N should
            # be R-N=[N+]=N
        ]

        # Add substructures and reactions (initially none)
        for i, rxn_datum in enumerate(rxn_data):
            rxn_data[i].append(Chem.MolFromSmarts(rxn_datum[0]))
            rxn_data[i].append(None)

        mol = smiles_to_mol(smi)
        # Add hydrogens (respects valence, so incomplete).
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol)

        while True:  # Keep going until all these issues have been resolved.
            current_rxn = None  # The reaction to perform.
            current_rxn_str = None

            for i, rxn_datum in enumerate(rxn_data):
                (
                    reactant_smarts,
                    product_smarts,
                    substruct_match_mol,
                    rxn_placeholder,
                ) = rxn_datum
                if mol.HasSubstructMatch(substruct_match_mol):
                    if rxn_placeholder is None:
                        current_rxn_str = reactant_smarts + ">>" + product_smarts
                        logger.trace("Defining reaction: {}", current_rxn_str)
                        current_rxn = AllChem.ReactionFromSmarts(current_rxn_str)
                        rxn_data[i][3] = current_rxn  # Update the placeholder.
                    else:
                        current_rxn = rxn_data[i][3]
                    break

            # Perform the reaction if necessary
            if current_rxn is None:  # No reaction left, so break out of while loop.
                break
            else:
                mol = current_rxn.RunReactants((mol,))[0][0]
                mol.UpdatePropertyCache(strict=False)  # Update valences

        # The mols have been altered from the reactions described above, we
        # need to resanitize them. Make sure aromatic rings are shown as such
        # This catches all RDKit Errors. without the catchError and
        # sanitizeOps the Chem.SanitizeMol can crash the program.
        sanitize_string = Chem.SanitizeMol(
            mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors=True
        )
        if sanitize_string.name == "SANITIZE_NONE":
            smiles_done = Chem.MolToSmiles(mol)
            logger.debug("Smiles after neutralization: {}", smiles_done)
            return smiles_done
        else:
            return None
