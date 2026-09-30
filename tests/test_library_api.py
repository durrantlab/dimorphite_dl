"""Checks the importable run() and run_with_mol_list() entry points."""

import sys
from pathlib import Path
from typing import Callable, Dict, List, Union

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


def test_run_returns_list(canonical_smiles: Callable[[str], str]) -> None:
    """Checks that run() passes main()'s result back to the caller."""

    output = dimorphite_dl.run(
        smiles="CCC(=O)O",
        min_ph=VERY_BASIC_PH,
        max_ph=VERY_BASIC_PH,
        return_as_list=True,
    )

    assert output is not None
    assert [line.split()[0] for line in output] == [canonical_smiles("CCC(=O)[O-]")]


def test_run_with_mol_list_returns_only_valid_mols(
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that a molecule with an unrepresentable protonated state does
    not put None in the returned list at the default pH range."""

    mols = dimorphite_dl.run_with_mol_list([Chem.MolFromSmiles("BrC1=CNC=C(C1=O)Br")])

    assert all(m is not None for m in mols)
    output = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
    assert output == [canonical_smiles("O=c1c(Br)c[nH]cc1Br")]


def test_run_with_mol_list_skips_none_entries(
    capsys: pytest.CaptureFixture[str],
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that a None from a failed Chem.MolFromSmiles is skipped with a
    warning, as bad SMILES are on the command line, instead of raising an
    opaque Boost error."""

    mols = dimorphite_dl.run_with_mol_list(
        [None, Chem.MolFromSmiles("CCCN")], min_ph=-1e7, max_ph=-1e7
    )

    output = [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
    assert output == [canonical_smiles("CCC[NH3+]")]
    assert "Skipping None entry" in capsys.readouterr().err


def test_protonate_does_not_modify_args() -> None:
    """Checks that a caller can reuse its args dict. Protonate used to add
    smiles_file to it, so a second call raised because both smiles and
    smiles_file were present."""

    args: Dict[str, Union[str, float]] = {
        "smiles": "CCCN",
        "min_ph": 7.0,
        "max_ph": 7.0,
    }
    original = dict(args)

    first = list(dimorphite_dl.Protonate(args))
    assert args == original

    second = list(dimorphite_dl.Protonate(args))
    assert second == first


@pytest.mark.parametrize(
    "params",
    [{}, {"smiles": "CCCN", "min_ph": 8.4, "max_ph": 6.4}],
    ids=["no_smiles", "inverted_ph_range"],
)
def test_invalid_arguments_leave_output_file_intact(
    params: Dict[str, Union[str, float]], tmp_path: Path
) -> None:
    """Checks that arguments are validated before the output file is opened
    for writing, so a bad call does not wipe earlier results."""

    output_file = tmp_path / "old.smi"
    output_file.write_text("previous results\n")

    with pytest.raises(Exception):
        dimorphite_dl.run(output_file=str(output_file), **params)

    assert output_file.read_text() == "previous results\n"
