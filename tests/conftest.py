"""Shared pytest setup for the Dimorphite-DL tests."""

import os
import sys
from typing import Callable

import pytest
from rdkit import Chem

# dimorphite_dl.py lives in the project root, which is not on sys.path when
# pytest imports tests from this directory.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _canonical_smiles(smiles: str) -> str:
    """Canonicalizes a SMILES string so expected values can be written in any
    valid form and still compare equal to Dimorphite-DL's output.

    Args:
        smiles: A valid SMILES string.

    Returns:
        RDKit's canonical isomeric SMILES for smiles.
    """

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


@pytest.fixture
def canonical_smiles() -> Callable[[str], str]:
    """Shares the canonicalizer across test modules without importing
    conftest directly.

    Returns:
        A function mapping a SMILES string to its canonical form.
    """

    return _canonical_smiles
