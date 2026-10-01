"""Checks how the command line and file reader handle input and output
files, and that stdout carries nothing but SMILES."""

import os
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

from dimorphite_dl import protonate_smiles

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# At the default pH range, the amine of CCCN is BOTH.
CCCN_STATES = ["CCCN", "CCC[NH3+]"]


def canonical(smiles: str) -> str:
    """Canonicalizes so expected values can be written in any valid form.

    Args:
        smiles: A valid SMILES string.

    Returns:
        RDKit's canonical isomeric SMILES for smiles.
    """

    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)


def run_cli(args: list[str], cwd: Path) -> "subprocess.CompletedProcess[str]":
    """Runs the command line in a fresh interpreter, because loguru binds its
    console sink when logging is configured, so this process cannot capture
    it reliably.

    Args:
        args: Command-line arguments after the program name.
        cwd: Working directory for the child process.

    Returns:
        The finished process, with stdout and stderr as text. A nonzero exit
        status does not raise.
    """

    env = dict(os.environ)
    env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    # Logging enabled from the environment would mix into every run.
    env.pop("DIMORPHITE_DL_LOG", None)
    return subprocess.run(
        [sys.executable, "-c", "from dimorphite_dl.cli import run_cli; run_cli()"]
        + args,
        cwd=str(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=False,
    )


def test_utf8_bom_does_not_drop_first_molecule(tmp_path: Path) -> None:
    """Checks that a BOM at the start of a file is ignored. Excel and Notepad
    add one, and it glued onto the first SMILES and made it unparseable."""

    path = tmp_path / "bom.smi"
    path.write_text("\ufeffCCCN\n", encoding="utf-8")

    output = protonate_smiles(str(path))

    assert sorted(canonical(s) for s in output) == sorted(
        canonical(s) for s in CCCN_STATES
    ), output


def test_log_lines_stay_out_of_stdout(tmp_path: Path) -> None:
    """Checks that --log_level sends log lines to stderr. They went to
    stdout, mixed in with the SMILES, so redirected output was corrupt."""

    result = run_cli(["--log_level", "debug", "CCCN"], tmp_path)

    assert result.returncode == 0, result.stderr
    assert result.stderr != ""

    # The point of the test is that log lines do not parse; keep RDKit from
    # printing an error for each one if they do appear.
    RDLogger.DisableLog("rdApp.*")
    try:
        lines = [line for line in result.stdout.splitlines() if line.strip()]
        parsed = [Chem.MolFromSmiles(line) for line in lines]
    finally:
        RDLogger.EnableLog("rdApp.*")

    assert all(mol is not None for mol in parsed), result.stdout
    assert sorted(canonical(line) for line in lines) == sorted(
        canonical(s) for s in CCCN_STATES
    )


def test_output_file_same_as_input_is_rejected(tmp_path: Path) -> None:
    """Checks that the command line refuses to write over its own input. The
    output file was opened for writing first, which truncated the input
    before it was read."""

    path = tmp_path / "molecules.smi"
    path.write_text("CCCN\n", encoding="utf-8")

    result = run_cli(["--output_file", str(path), str(path)], tmp_path)

    assert result.returncode != 0
    assert "is the input file" in result.stderr, result.stderr
    assert path.read_text(encoding="utf-8") == "CCCN\n"


def test_invalid_arguments_leave_output_file_intact(tmp_path: Path) -> None:
    """Checks that arguments rejected during protonation do not truncate an
    existing output file, which was opened before anything was checked."""

    path = tmp_path / "existing.smi"
    path.write_text("keep me\n", encoding="utf-8")

    result = run_cli(
        ["--ph_min", "9", "--ph_max", "5", "--output_file", str(path), "CCCN"],
        tmp_path,
    )

    assert result.returncode != 0
    assert path.read_text(encoding="utf-8") == "keep me\n"


def test_output_file_is_written(tmp_path: Path) -> None:
    """Checks that --output_file still receives every state. The file is now
    written in a with block rather than left for the interpreter to flush."""

    path = tmp_path / "out.smi"

    result = run_cli(["--output_file", str(path), "CCCN"], tmp_path)

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    lines = path.read_text(encoding="utf-8").splitlines()
    assert sorted(canonical(line) for line in lines) == sorted(
        canonical(s) for s in CCCN_STATES
    )
