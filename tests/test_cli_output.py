"""Checks that the command line writes nothing but SMILES to stdout, so its
output can be redirected straight into a .smi file."""

import os
import subprocess
import sys
from pathlib import Path
from typing import Callable, List

from rdkit import Chem

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(PROJECT_ROOT, "dimorphite_dl.py")


def run_python(args: List[str], cwd: str) -> str:
    """Runs a fresh interpreter, because the banner prints at import time and
    this test process has already imported the module.

    Args:
        args: Arguments passed to the Python interpreter.
        cwd: Working directory for the child process.

    Returns:
        Everything the child wrote to stdout.
    """

    result = subprocess.run(
        [sys.executable] + args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )
    return result.stdout


def test_command_line_stdout_is_only_smiles(
    tmp_path: Path, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that the banner and PARAMETERS block stay out of stdout."""

    stdout = run_python(
        [SCRIPT, "--smiles", "CCCN", "--min_ph=-10000000", "--max_ph=-10000000"],
        str(tmp_path),
    )

    lines = [line for line in stdout.splitlines() if line.strip() != ""]
    assert len(lines) == 1, stdout
    assert Chem.MolFromSmiles(lines[0].split()[0]) is not None, stdout
    assert lines[0].split()[0] == canonical_smiles("CCC[NH3+]"), stdout


def test_import_writes_nothing_to_stdout() -> None:
    """Checks that importing the module leaves the host's stdout clean."""

    assert run_python(["-c", "import dimorphite_dl"], PROJECT_ROOT) == ""
