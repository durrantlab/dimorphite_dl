from typing import Any

import os
import sys
from io import StringIO

from rdkit import Chem
from rdkit.Chem import AllChem


def convert_smiles_str_to_mol(smiles_str):
    """Given a SMILES string, check that it is actually a string and not a
    None. Then try to convert it to an RDKit Mol Object.

    Args:
        smiles_str: The SMILES string.

    Returns:
        A rdkit.Chem.rdchem.Mol object, or None if it is the wrong type or
            if it fails to convert to a Mol Obj
    """

    # Check that there are no type errors, ie Nones or non-string A
    # non-string type will cause RDKit to hard crash
    if smiles_str is None or type(smiles_str) is not str:
        return None

    # Try to fix azides here. They are just tricky to deal with.
    smiles_str = smiles_str.replace("N=N=N", "N=[N+]=N")
    smiles_str = smiles_str.replace("NN#N", "N=[N+]=N")

    # Now convert to a mol object. Note the trick that is necessary to
    # capture RDKit error/warning messages. See
    # https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable
    stderr_fileno = sys.stderr.fileno()
    stderr_save = os.dup(stderr_fileno)
    stderr_pipe = os.pipe()
    os.dup2(stderr_pipe[1], stderr_fileno)
    os.close(stderr_pipe[1])

    mol = Chem.MolFromSmiles(smiles_str)

    os.close(stderr_fileno)
    os.close(stderr_pipe[0])
    os.dup2(stderr_save, stderr_fileno)
    os.close(stderr_save)

    # Check that there are None type errors Chem.MolFromSmiles has
    # sanitize on which means if there is even a small error in the SMILES
    # (kekulize, nitrogen charge...) then mol=None. ie.
    # Chem.MolFromSmiles("C[N]=[N]=[N]") = None this is an example of an
    # nitrogen charge error. It is cased in a try statement to be overly
    # cautious.
    return None if mol is None else mol


def neutralize_mol(mol):
    """All molecules should be neuralized to the extent possible. The user
    should not be allowed to specify the valence of the atoms in most cases.

    Args:
        mol: The rdkit Mol objet to be neutralized.

    Returns:
        The neutralized Mol object.
    """

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

    # Add hydrogens (respects valence, so incomplete).
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol)

    while True:  # Keep going until all these issues have been resolved.
        current_rxn = None  # The reaction to perform.
        current_rxn_str = None

        for i, rxn_datum in enumerate(rxn_data):
            (reactant_smarts, product_smarts, substruct_match_mol, rxn_placeholder) = (
                rxn_datum
            )
            if mol.HasSubstructMatch(substruct_match_mol):
                if rxn_placeholder is None:
                    current_rxn_str = reactant_smarts + ">>" + product_smarts
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

    return mol if sanitize_string.name == "SANITIZE_NONE" else None


class LoadSMIFile:
    """A generator class for loading in the SMILES strings from a file, one at
    a time."""

    def __init__(self, filename: StringIO | str, kwargs: dict[str, Any]) -> None:
        """Initializes this class.

        Args:
            filename: The filename or file object (i.e., StringIO).
        """

        self.kwargs = kwargs

        if type(filename) is str:
            # It's a filename
            self.f = open(filename, "r")
        else:
            # It's a file object (i.e., StringIO)
            self.f = filename

    def __iter__(self):
        """Returns this generator object."""
        return self

    def __next__(self) -> dict[str, str]:
        """A dict, where the "smiles" key contains the canonical SMILES
        string and the "data" key contains the remaining information
        (e.g., the molecule name).
        """
        return self.next()

    def next(self) -> dict[str, Any]:
        """Get the data associated with the next line.

        Returns:
            A dict, where the "smiles" key contains the canonical SMILES
                 string and the "data" key contains the remaining information
                 (e.g., the molecule name).
        """

        line = self.f.readline()

        if line == "":
            # EOF
            self.f.close()
            raise StopIteration()

        # Divide line into smi and data
        splits = line.split()
        if len(splits) != 0:
            # Generate mol object
            smiles_str = splits[0]

            # Convert from SMILES string to RDKIT Mol. This series of tests is
            # to make sure the SMILES string is properly formed and to get it
            # into a canonical form. Filter if failed.
            mol = convert_smiles_str_to_mol(smiles_str)
            if mol is None:
                if "silent" in self.kwargs and not self.kwargs["silent"]:
                    print("WARNING: Skipping poorly formed SMILES string: " + line)
                return self.next()

            # Handle nuetralizing the molecules. Filter if failed.
            mol = neutralize_mol(mol)
            if mol is None:
                if "silent" in self.kwargs and not self.kwargs["silent"]:
                    print("WARNING: Skipping poorly formed SMILES string: " + line)
                return self.next()

            # Remove the hydrogens.
            try:
                mol = Chem.RemoveHs(mol)
            except:
                if "silent" in self.kwargs and not self.kwargs["silent"]:
                    print("WARNING: Skipping poorly formed SMILES string: " + line)
                return self.next()

            if mol is None:
                if "silent" in self.kwargs and not self.kwargs["silent"]:
                    print("WARNING: Skipping poorly formed SMILES string: " + line)
                return self.next()

            # Regenerate the smiles string (to standardize).
            new_mol_string = Chem.MolToSmiles(mol, isomericSmiles=True)

            return {"smiles": new_mol_string, "data": splits[1:]}
        else:
            # Blank line? Go to next one.
            return self.next()
