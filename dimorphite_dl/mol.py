import os
import sys

from loguru import logger
from rdkit import Chem


def smiles_to_mol(smiles_str: str) -> Chem.Mol | None:
    """Convert a SMILES string to an RDKit Mol object with detailed error logging.

    Args:
        smiles_str: The SMILES string to convert.

    Returns:
        RDKit Mol object if successful, None if conversion fails.

    Raises:
        TypeError: If input is not a string.
    """
    if not smiles_str.strip():
        logger.error("SMILES input is empty or whitespace only")
        return None

    # Clean up the input
    smiles_str = smiles_str.strip()
    logger.debug("Converting SMILES to RDKit mol: '{}'", smiles_str)

    # Handle common azide patterns
    original_smiles = smiles_str
    if "N=N=N" in smiles_str or "NN#N" in smiles_str:
        logger.info("Attempting to fix azide patterns in: '{}'", smiles_str)
        smiles_str = smiles_str.replace("N=N=N", "N=[N+]=N")
        smiles_str = smiles_str.replace("NN#N", "N=[N+]=N")
        if smiles_str != original_smiles:
            logger.info("Modified SMILES: '{}' -> '{}'", original_smiles, smiles_str)

    # Capture RDKit stderr to get detailed error messages
    rdkit_errors = _capture_rdkit_errors(smiles_str)

    if rdkit_errors["mol"] is None:
        error_msg = rdkit_errors["stderr_content"].strip()
        if error_msg:
            logger.warning(
                "RDKit failed to parse SMILES '{}'. RDKit error: {}",
                smiles_str,
                error_msg,
            )
        else:
            logger.warning(
                "RDKit failed to parse SMILES '{}' (no specific error message)",
                smiles_str,
            )
        return None

    logger.debug("Successfully converted SMILES to RDKit mol object")
    return rdkit_errors["mol"]


def _capture_rdkit_errors(smiles_str: str) -> dict:
    """Capture RDKit stderr output and return mol object with error messages.

    Args:
        smiles_str: SMILES string to convert

    Returns:
        Dictionary with 'mol' (RDKit Mol or None) and 'stderr_content' (string)
    """
    # Set up stderr capture
    stderr_fileno = sys.stderr.fileno()
    stderr_save = os.dup(stderr_fileno)
    stderr_pipe = os.pipe()

    try:
        # Redirect stderr to pipe
        os.dup2(stderr_pipe[1], stderr_fileno)
        os.close(stderr_pipe[1])

        # Convert SMILES to mol (this may write to stderr)
        mol = Chem.MolFromSmiles(smiles_str)

        # Read captured stderr
        os.close(stderr_fileno)
        stderr_content = os.read(stderr_pipe[0], 1024).decode("utf-8", errors="ignore")

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
