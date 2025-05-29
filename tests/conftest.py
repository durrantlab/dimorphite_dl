import os

import pytest
from rdkit import Chem

from dimorphite_dl import enable_logging
from dimorphite_dl.io import SMILESProcessor

TEST_DIR = os.path.dirname(__file__)


def compare_smiles(smiles1, smiles2):
    detected_can = Chem.MolToSmiles(Chem.MolFromSmiles(smiles1), isomericSmiles=True)
    expected_can = Chem.MolToSmiles(Chem.MolFromSmiles(smiles2), isomericSmiles=True)
    assert detected_can == expected_can, f"got {smiles1}, expected {smiles2}"


def compare_smarts(smarts1, smarts2):
    detected_can = Chem.MolToSmarts(Chem.MolFromSmarts(smarts1))
    expected_can = Chem.MolToSmarts(Chem.MolFromSmarts(smarts2))
    assert detected_can == expected_can, f"got {smarts1}, expected {smarts2}"


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)


# Pytest fixtures for reusable test data
@pytest.fixture
def sample_smiles_list():
    """Fixture providing a sample list of SMILES strings."""
    return ["CCO", "CCC", "c1ccccc1", "CC(C)C", "CCN"]


@pytest.fixture
def sample_smiles_file(tmp_path):
    """Fixture providing a temporary SMILES file."""
    content = "CCO ethanol\nCCC propane\nc1ccccc1 benzene\n"
    file_path = tmp_path / "test_molecules.smi"
    file_path.write_text(content)
    return str(file_path)


@pytest.fixture
def processor_no_validation():
    """Fixture providing a processor with validation disabled."""
    return SMILESProcessor(validate_smiles=False)
