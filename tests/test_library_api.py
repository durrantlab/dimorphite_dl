"""Checks the importable run() and run_with_mol_list() entry points."""

import sys
from pathlib import Path
from typing import Callable, List

import pytest
from rdkit import Chem

import dimorphite_dl

VERY_BASIC_PH = 10000000.0

# Flags a host process might have in sys.argv: a Jupyter kernel's connection
# file, and a host script's own option that shares a name with ours.
HOST_FLAGS = [["-f", "kernel.json"], ["--output_file", "stray.smi"]]
HOST_FLAG_IDS = ["jupyter_kernel", "host_output_file"]


@pytest.mark.parametrize("host_flags", HOST_FLAGS, ids=HOST_FLAG_IDS)
def test_run_with_mol_list_ignores_sys_argv(
    host_flags: List[str],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that library calls return their result instead of parsing the
    host's command line or writing a file it names."""

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["host_script.py"] + host_flags)

    mols = dimorphite_dl.run_with_mol_list(
        [Chem.MolFromSmiles("CCC(=O)O")], min_ph=VERY_BASIC_PH, max_ph=VERY_BASIC_PH
    )

    output = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
    assert output == [canonical_smiles("CCC(=O)[O-]")]
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("host_flags", HOST_FLAGS, ids=HOST_FLAG_IDS)
def test_run_ignores_sys_argv(
    host_flags: List[str],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that run() prints its result instead of parsing the host's
    command line or writing a file it names."""

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["host_script.py"] + host_flags)

    dimorphite_dl.run(smiles="CCC(=O)O", min_ph=VERY_BASIC_PH, max_ph=VERY_BASIC_PH)

    assert canonical_smiles("CCC(=O)[O-]") in capsys.readouterr().out.split()
    assert list(tmp_path.iterdir()) == []
