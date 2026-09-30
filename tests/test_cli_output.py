"""Checks that the command line writes nothing but SMILES to stdout, so its
output can be redirected straight into a .smi file."""

import os
import subprocess
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional

from rdkit import Chem

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(PROJECT_ROOT, "dimorphite_dl.py")


def run_python(
    args: List[str],
    cwd: str,
    env: Optional[Dict[str, str]] = None,
    check: bool = True,
) -> "subprocess.CompletedProcess[str]":
    """Runs a fresh interpreter, because import-time behavior cannot be
    observed from this test process, which has already imported the module.

    Args:
        args: Arguments passed to the Python interpreter.
        cwd: Working directory for the child process.
        env: Environment for the child process. Defaults to this process's.
        check: Whether a nonzero exit status raises CalledProcessError.

    Returns:
        The finished process, with its stdout and stderr as text.
    """

    return subprocess.run(
        [sys.executable] + args,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=check,
    )


def test_command_line_stdout_is_only_smiles(
    tmp_path: Path, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that the banner and PARAMETERS block stay out of stdout."""

    result = run_python(
        [SCRIPT, "--smiles", "CCCN", "--min_ph=-10000000", "--max_ph=-10000000"],
        str(tmp_path),
    )
    stdout = result.stdout

    lines = [line for line in stdout.splitlines() if line.strip() != ""]
    assert len(lines) == 1, stdout
    assert Chem.MolFromSmiles(lines[0].split()[0]) is not None, stdout
    assert lines[0].split()[0] == canonical_smiles("CCC[NH3+]"), stdout

    # The command line still shows the help hint and citation.
    assert "please cite" in result.stderr, result.stderr


def test_argument_error_keeps_stdout_empty(tmp_path: Path) -> None:
    """Checks that a bad argument writes its help text and error to stderr.
    They went to stdout, so `> out.smi` put them in the file, where a
    downstream tool would read them as SMILES."""

    result = run_python(
        [SCRIPT, "--smiles", "CCCN", "--min_ph", "abc"], str(tmp_path), check=False
    )

    assert result.returncode != 0
    assert result.stdout == ""
    assert "ERROR:" in result.stderr, result.stderr
    assert "examples:" in result.stderr, result.stderr


def test_output_order_is_independent_of_hash_seed(tmp_path: Path) -> None:
    """Checks that variant order is stable across interpreters, because
    callers that keep only the first N variants would otherwise get different
    molecules on identical reruns. String hashing is randomized per process,
    so any set-based deduplication would reorder the output."""

    # Two BOTH sites give four variants, enough for a reordering to show.
    args = [
        SCRIPT,
        "--smiles",
        "NCCCC(=O)O",
        "--min_ph=7.0",
        "--max_ph=7.0",
        "--pka_precision=1000000.0",
    ]

    outputs = []
    for seed in ["1", "2", "3"]:
        env = dict(os.environ, PYTHONHASHSEED=seed)
        outputs.append(run_python(args, str(tmp_path), env).stdout)

    assert len(outputs[0].splitlines()) == 4, outputs[0]
    assert outputs[1] == outputs[0]
    assert outputs[2] == outputs[0]


def test_import_writes_nothing() -> None:
    """Checks that importing the module prints nothing to either stream, so
    it can be embedded in other programs."""

    result = run_python(["-c", "import dimorphite_dl"], PROJECT_ROOT)

    assert result.stdout == ""
    assert result.stderr == ""


def test_missing_rdkit_raises_import_error() -> None:
    """Checks that a missing RDKit raises ImportError, which callers can
    catch specifically, and prints nothing to stdout."""

    # Setting a sys.modules entry to None makes any import of it fail.
    script = (
        "import sys\n"
        "sys.modules['rdkit'] = None\n"
        "try:\n"
        "    import dimorphite_dl\n"
        "except ImportError as e:\n"
        "    print('ImportError: ' + str(e))\n"
    )

    result = run_python(["-c", script], PROJECT_ROOT)

    assert result.stdout == (
        "ImportError: Dimorphite-DL requires RDKit. See https://www.rdkit.org/\n"
    )
