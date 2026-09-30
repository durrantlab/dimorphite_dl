"""Checks how input SMILES are read, skipped, and validated before any
protonation happens."""

import warnings
from io import StringIO
from pathlib import Path

import pytest
from rdkit import Chem

import dimorphite_dl


def test_long_run_of_skipped_lines(capsys: pytest.CaptureFixture[str]) -> None:
    """Checks that thousands of blank and unparseable lines in a row are
    skipped without exhausting the stack."""

    text = "\n" * 3000 + "not_a_smiles\n" * 3000 + "CCCN name\n"

    records = list(dimorphite_dl.LoadSMIFile(StringIO(text)))

    assert records == [{"smiles": "CCCN", "data": ["name"]}]
    assert "Skipping poorly formed SMILES string" in capsys.readouterr().err


def test_rdkit_error_skips_the_line(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Checks that an RDKit failure on one molecule skips only that line."""

    def failing_remove_hs(mol: Chem.Mol) -> Chem.Mol:
        """Raises what RDKit raises for a molecule it cannot sanitize.

        Args:
            mol: Ignored.

        Returns:
            Never returns.
        """

        raise ValueError("Sanitization error")

    monkeypatch.setattr(dimorphite_dl.Chem, "RemoveHs", failing_remove_hs)

    assert list(dimorphite_dl.LoadSMIFile(StringIO("CCCN\n"))) == []
    assert "Skipping poorly formed SMILES string" in capsys.readouterr().err


def test_programming_error_is_not_reported_as_bad_smiles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Checks that a bug inside the pipeline propagates instead of being
    turned into a skipped input."""

    def broken_remove_hs(mol: Chem.Mol) -> Chem.Mol:
        """Stands in for a coding error, not an RDKit failure.

        Args:
            mol: Ignored.

        Returns:
            Never returns.
        """

        raise AttributeError("simulated bug")

    monkeypatch.setattr(dimorphite_dl.Chem, "RemoveHs", broken_remove_hs)

    with pytest.raises(AttributeError):
        list(dimorphite_dl.LoadSMIFile(StringIO("CCCN\n")))


def test_path_object_is_read_as_a_filename(tmp_path: Path) -> None:
    """Checks that a pathlib.Path is opened like a str path. Only str was
    treated as a filename, so a Path was used as a file object and crashed on
    readline()."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCCN\n")

    assert list(dimorphite_dl.LoadSMIFile(smi_file)) == [{"smiles": "CCCN", "data": []}]


def test_smiles_and_smiles_file_together_are_rejected(tmp_path: Path) -> None:
    """Checks that giving both inputs raises instead of silently ignoring the
    file."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCC(=O)O\n")

    with pytest.raises(ValueError):
        dimorphite_dl.Protonate({"smiles": "CCCN", "smiles_file": str(smi_file)})


def test_module_source_has_no_syntax_warnings() -> None:
    """Checks for comparisons like `is not ""`, which only work because of
    string interning and raise SyntaxWarning on Python 3.8+."""

    with open(dimorphite_dl.__file__) as source_file:
        source = source_file.read()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compile(source, dimorphite_dl.__file__, "exec")
