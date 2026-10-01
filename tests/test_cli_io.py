"""Checks how the command line and file reader handle input and output
files, and that stdout carries nothing but SMILES."""

import os
import re
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest
from loguru import logger
from rdkit import Chem, RDLogger

from dimorphite_dl import enable_logging, protonate_smiles
from dimorphite_dl.io import SMILESStreamError

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


def run_cli(
    args: list[str], cwd: Path, env_overrides: dict[str, str] | None = None
) -> "subprocess.CompletedProcess[str]":
    """Runs the command line in a fresh interpreter, because loguru binds its
    console sink when logging is configured, so this process cannot capture
    it reliably.

    Args:
        args: Command-line arguments after the program name.
        cwd: Working directory for the child process.
        env_overrides: Extra environment variables, such as HOME for tests
            that pass a "~" path.

    Returns:
        The finished process, with stdout and stderr as text. A nonzero exit
        status does not raise.
    """

    env = dict(os.environ)
    env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    # Logging enabled from the environment would mix into every run.
    env.pop("DIMORPHITE_DL_LOG", None)
    env.update(env_overrides or {})
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


def test_identifiers_are_kept(tmp_path: Path) -> None:
    """Checks that names from the input file reach the output. They were
    dropped, so variants could not be traced back to their input."""

    path = tmp_path / "molecules.smi"
    path.write_text("CCCN amine\nCCO ethanol\n", encoding="utf-8")

    result = run_cli([str(path)], tmp_path)

    assert result.returncode == 0, result.stderr
    pairs = [line.split(",") for line in result.stdout.splitlines()]
    assert sorted(canonical(smi) for smi, name in pairs if name == "amine") == sorted(
        canonical(s) for s in CCCN_STATES
    ), result.stdout
    assert [canonical(smi) for smi, name in pairs if name == "ethanol"] == [
        canonical("CCO")
    ], result.stdout


def test_missing_input_file_is_an_error(tmp_path: Path) -> None:
    """Checks that a missing input file fails with a message and leaves an
    existing output file alone. The error was swallowed and the run exited 0
    with no output."""

    out = tmp_path / "existing.smi"
    out.write_text("keep me\n", encoding="utf-8")

    result = run_cli(
        ["--output_file", str(out), str(tmp_path / "missing.smi")], tmp_path
    )

    assert result.returncode != 0
    assert "File not found" in result.stderr, result.stderr
    assert out.read_text(encoding="utf-8") == "keep me\n"


def test_output_is_written_as_it_is_produced(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Checks that the command line writes each result as it is produced.
    Results were collected into a list first, so memory grew with the
    library and a failure late in the run left no output at all."""

    from dimorphite_dl import cli

    def failing_protonator(**kwargs: object) -> Iterator[str]:
        yield "CCCN"
        raise SMILESStreamError("input ended unexpectedly")

    out = tmp_path / "out.smi"
    monkeypatch.setattr(cli, "Protonate", failing_protonator)
    monkeypatch.setattr(
        sys, "argv", ["dimorphite_dl", "--output_file", str(out), "CCCN"]
    )

    with pytest.raises(SystemExit) as excinfo:
        cli.run_cli()

    assert excinfo.value.code == 1
    assert out.read_text(encoding="utf-8") == "CCCN\n"


def test_env_flag_values_do_not_break_import() -> None:
    """Checks that common spellings of the logging environment variables are
    accepted. literal_eval and int() raised on "true" and "INFO", so the
    package could not be imported."""

    env = dict(os.environ)
    env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env["DIMORPHITE_DL_LOG"] = "true"
    env["DIMORPHITE_DL_LOG_LEVEL"] = "INFO"
    env["DIMORPHITE_DL_STDOUT"] = "yes"
    result = subprocess.run(
        [sys.executable, "-c", "import dimorphite_dl"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_output_file_same_as_tilde_input_is_rejected(tmp_path: Path) -> None:
    """Checks that a quoted "~" input is still recognized as the output file.
    The reader expanded "~" but the guard did not, so the guard was skipped
    and the input was truncated while it was being read."""

    path = tmp_path / "molecules.smi"
    path.write_text("CCCN\nCCO\n", encoding="utf-8")
    home = {"HOME": str(tmp_path), "USERPROFILE": str(tmp_path)}

    result = run_cli(["--output_file", str(path), "~/molecules.smi"], tmp_path, home)

    assert result.returncode != 0
    assert "is the input file" in result.stderr, result.stderr
    assert path.read_text(encoding="utf-8") == "CCCN\nCCO\n"


def test_tilde_output_file_is_expanded(tmp_path: Path) -> None:
    """Checks that a quoted "~" output path lands in the home directory. It
    was opened literally, which failed because "./~/" does not exist."""

    home_dir = tmp_path / "home"
    home_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    home = {"HOME": str(home_dir), "USERPROFILE": str(home_dir)}

    result = run_cli(["--output_file", "~/out.smi", "CCCN"], work_dir, home)

    assert result.returncode == 0, result.stderr
    assert not (work_dir / "~").exists()
    lines = (home_dir / "out.smi").read_text(encoding="utf-8").splitlines()
    assert sorted(canonical(line) for line in lines) == sorted(
        canonical(s) for s in CCCN_STATES
    )


def test_log_file_has_no_color_codes(tmp_path: Path) -> None:
    """Checks that log files are plain text. colorize was passed to the file
    sink too, so ANSI escape codes were written into the file."""

    path = tmp_path / "run.log"
    enable_logging(20, stdout_set=False, file_path=str(path))
    try:
        logger.info("probe message")
    finally:
        # Matches the session fixture in conftest.py; also closes the file.
        enable_logging(0)

    text = path.read_text(encoding="utf-8")
    assert "probe message" in text
    assert "\x1b[" not in text, text


def test_enable_logging_keeps_host_sinks() -> None:
    """Checks that enabling logging leaves the host application's loguru
    sinks alone. logger.configure replaced every handler, so the host's sink
    stopped receiving messages and removing it raised ValueError."""

    messages: list[str] = []
    host_id = logger.add(messages.append, format="{message}")
    try:
        enable_logging(20, stdout_set=False)
        logger.info("host probe")
    finally:
        # Matches the session fixture in conftest.py.
        enable_logging(0)
    logger.remove(host_id)

    assert any("host probe" in message for message in messages), messages


@pytest.mark.parametrize(
    ("args", "env_overrides"),
    [
        (["--log_level", "info", "CCCN"], {}),
        (["CCCN"], {"DIMORPHITE_DL_LOG": "1", "DIMORPHITE_DL_LOG_LEVEL": "INFO"}),
    ],
)
def test_cli_logs_each_line_once(
    tmp_path: Path, args: list[str], env_overrides: dict[str, str]
) -> None:
    """Checks that the command line does not also log through loguru's
    default sink, which enable_logging no longer removes. Its lines start
    with a full date, unlike LOG_FORMAT."""

    result = run_cli(args, tmp_path, env_overrides)

    assert result.returncode == 0, result.stderr
    lines = [re.sub(r"\x1b\[[0-9;]*m", "", line) for line in result.stderr.splitlines()]
    assert lines, result.stderr
    default_format = [line for line in lines if re.match(r"\d{4}-\d{2}-\d{2} ", line)]
    assert default_format == [], result.stderr


def test_unreadable_input_is_an_error(tmp_path: Path) -> None:
    """Checks that input with no readable molecules fails with a message and
    leaves an existing output file alone. A mistyped path such as "input" is
    read as a SMILES string, rejected, and the run exited 0 with no output."""

    out = tmp_path / "existing.smi"
    out.write_text("keep me\n", encoding="utf-8")

    result = run_cli(["--output_file", str(out), "input"], tmp_path)

    assert result.returncode == 1, result.stderr
    assert "no readable molecules" in result.stderr, result.stderr
    assert out.read_text(encoding="utf-8") == "keep me\n"


def test_skipped_lines_are_reported(tmp_path: Path) -> None:
    """Checks that an unreadable line is reported on stderr while the rest
    of the file is still processed. It was dropped without a word, so the
    output silently had fewer molecules than the input."""

    path = tmp_path / "molecules.smi"
    path.write_text("CCO\nnot_a_smiles\n", encoding="utf-8")

    result = run_cli([str(path)], tmp_path)

    assert result.returncode == 0, result.stderr
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    assert [canonical(line) for line in lines] == [canonical("CCO")], result.stdout
    assert "skipped 1 unreadable input line(s)" in result.stderr, result.stderr


def test_clean_run_writes_nothing_to_stderr(tmp_path: Path) -> None:
    """Checks that the new warnings stay quiet when every input is
    protonated, so stderr remains usable as a failure signal."""

    result = run_cli(["CCCN"], tmp_path)

    assert result.returncode == 0, result.stderr
    assert result.stderr == ""


def test_unprotonated_fallback_is_reported(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Checks that a molecule written out unchanged after a protonation error
    is reported. The fallback line looks exactly like a real result and was
    only logged, which is off by default."""

    from dimorphite_dl import cli
    from dimorphite_dl.protonate.run import Protonate

    def failing_variants(*args: object, **kwargs: object) -> list[Chem.Mol]:
        raise RuntimeError("simulated failure")

    monkeypatch.setattr(Protonate, "_generate_protonated_variants", failing_variants)
    monkeypatch.setattr(sys, "argv", ["dimorphite_dl", "CCCN"])

    cli.run_cli()

    captured = capsys.readouterr()
    assert [canonical(line) for line in captured.out.splitlines()] == [
        canonical("CCCN")
    ], captured.out
    assert "1 molecule(s) could not be protonated" in captured.err, captured.err


def test_out_of_range_arguments_raise_value_error() -> None:
    """Checks that bad ranges raise ValueError. They were checked with
    assert, which python -O strips, so an inverted pH range ran anyway."""

    with pytest.raises(ValueError, match="ph_min"):
        protonate_smiles("CCCN", ph_min=9.0, ph_max=5.0)
    with pytest.raises(ValueError, match="precision"):
        protonate_smiles("CCCN", precision=-1.0)
    with pytest.raises(ValueError, match="max_variants"):
        protonate_smiles("CCCN", max_variants=0)


def test_invalid_arguments_are_a_usage_error(tmp_path: Path) -> None:
    """Checks that an inverted pH range on the command line gives a usage
    message rather than a traceback."""

    result = run_cli(["--ph_min", "9", "--ph_max", "5", "CCCN"], tmp_path)

    assert result.returncode == 2, result.stderr
    assert "must be less than or equal to ph_max" in result.stderr, result.stderr
    assert "Traceback" not in result.stderr, result.stderr
