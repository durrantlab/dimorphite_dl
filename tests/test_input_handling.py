"""Checks how input SMILES are read, skipped, and validated before any
protonation happens."""

import inspect
import warnings
from io import StringIO
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

import dimorphite_dl


def test_long_run_of_skipped_lines(capsys: pytest.CaptureFixture[str]) -> None:
    """Checks that thousands of blank and unparseable lines in a row are
    skipped without exhausting the stack."""

    text = "\n" * 3000 + "not_a_smiles\n" * 3000 + "CCCN name\n"

    records = list(dimorphite_dl.LoadSMIFile(StringIO(text)))

    assert records == [{"smiles": "CCCN", "data": ["name"]}]
    assert "RDKit could not parse it" in capsys.readouterr().err


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
    assert "removing hydrogens failed: Sanitization error" in capsys.readouterr().err


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


def test_input_file_closed_when_protonation_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Checks that main() closes the input file when an exception stops it
    partway through. The file was closed only at EOF, so each failed call
    from a long-running library caller leaked a handle."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCCN\nCCC(=O)O\n")

    loaders: List[dimorphite_dl.LoadSMIFile] = []
    original_init = dimorphite_dl.LoadSMIFile.__init__

    def recording_init(self: dimorphite_dl.LoadSMIFile, filename: str) -> None:
        """Keeps each loader so the test can inspect its file afterward.

        Args:
            self: The loader being constructed.
            filename: Passed through unchanged.
        """

        original_init(self, filename)
        loaders.append(self)

    def failing_protonate_sites(*args: object) -> None:
        """Stands in for an error raised partway through the input.

        Args:
            *args: Ignored.
        """

        raise RuntimeError("simulated failure")

    monkeypatch.setattr(dimorphite_dl.LoadSMIFile, "__init__", recording_init)
    monkeypatch.setattr(
        dimorphite_dl.ProtSubstructFuncs,
        "protonate_sites",
        staticmethod(failing_protonate_sites),
    )

    with pytest.raises(RuntimeError, match="simulated failure"):
        dimorphite_dl.run(smiles_file=str(smi_file), return_as_list=True)

    assert len(loaders) == 1
    assert loaders[0].f.closed


def test_input_file_closed_when_substructure_loading_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Checks that Protonate closes the input file when loading the SMARTS
    substructures fails. The file is opened before loading, and main() cannot
    close it because the constructor never returns."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCCN\n")

    loaders: List[dimorphite_dl.LoadSMIFile] = []
    original_init = dimorphite_dl.LoadSMIFile.__init__

    def recording_init(self: dimorphite_dl.LoadSMIFile, filename: str) -> None:
        """Keeps each loader so the test can inspect its file afterward.

        Args:
            self: The loader being constructed.
            filename: Passed through unchanged.
        """

        original_init(self, filename)
        loaders.append(self)

    def failing_loader(*args: object) -> None:
        """Stands in for a missing or malformed SMARTS file.

        Args:
            *args: Ignored.
        """

        raise RuntimeError("simulated failure")

    monkeypatch.setattr(dimorphite_dl.LoadSMIFile, "__init__", recording_init)
    monkeypatch.setattr(
        dimorphite_dl.ProtSubstructFuncs,
        "load_protonation_substructs_calc_state_for_ph",
        staticmethod(failing_loader),
    )

    with pytest.raises(RuntimeError, match="simulated failure"):
        dimorphite_dl.run(smiles_file=str(smi_file), return_as_list=True)

    assert len(loaders) == 1
    assert loaders[0].f.closed


def test_load_smi_file_closes_on_early_exit(tmp_path: Path) -> None:
    """Checks that a caller who stops iterating before EOF can still close
    the file, and that the later EOF close does not fail."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCCN\nCCC(=O)O\n")

    with dimorphite_dl.LoadSMIFile(smi_file) as loader:
        assert loader.next()["smiles"] == "CCCN"
    assert loader.f.closed

    loader.close()


def test_smiles_and_smiles_file_together_are_rejected(tmp_path: Path) -> None:
    """Checks that giving both inputs raises instead of silently ignoring the
    file."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_text("CCC(=O)O\n")

    with pytest.raises(ValueError):
        dimorphite_dl.Protonate({"smiles": "CCCN", "smiles_file": str(smi_file)})


def test_missing_smiles_error_goes_to_stderr(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Checks that the "No SMILES" error stays off stdout, which carries only
    protonated SMILES."""

    with pytest.raises(Exception, match="No SMILES"):
        dimorphite_dl.Protonate({})

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "No SMILES" in captured.err


def test_module_source_has_no_syntax_warnings() -> None:
    """Checks for comparisons like `is not ""`, which only work because of
    string interning and raise SyntaxWarning on Python 3.8+."""

    with open(dimorphite_dl.__file__) as source_file:
        source = source_file.read()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        compile(source, dimorphite_dl.__file__, "exec")


def test_utf8_bom_does_not_drop_first_molecule(tmp_path: Path) -> None:
    """Checks that the BOM Excel and Notepad write is stripped. It stuck to
    the first SMILES, which was then skipped as poorly formed."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_bytes(b"\xef\xbb\xbfCCCN name\n")

    assert list(dimorphite_dl.LoadSMIFile(smi_file)) == [
        {"smiles": "CCCN", "data": ["name"]}
    ]


def test_non_ascii_names_round_trip_as_utf8(tmp_path: Path) -> None:
    """Checks that input and output files are UTF-8 regardless of the locale.
    With the locale default (cp1252 on Windows), non-ASCII names were garbled
    or raised UnicodeDecodeError mid-run."""

    smi_file = tmp_path / "input.smi"
    smi_file.write_bytes("CCCC caf\u00e9_\u03b1\n".encode("utf-8"))
    output_file = tmp_path / "output.smi"

    dimorphite_dl.run(smiles_file=str(smi_file), output_file=str(output_file))

    assert output_file.read_bytes().decode("utf-8").split() == [
        "CCCC",
        "caf\u00e9_\u03b1",
    ]


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(=O)O[2H]", "CC(=O)O"),
        ("CC(=O)O[3H]", "CC(=O)O"),
        ("[2H]C([2H])([2H])C(=O)O[2H]", "[2H]C([2H])([2H])C(=O)O"),
    ],
    ids=["deuterium", "tritium", "carbon_label_kept"],
)
def test_exchangeable_h_isotopes_become_plain_h(
    smiles: str, expected: str, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that D and T on heteroatoms are loaded as ordinary hydrogens
    while labels on carbon survive. RemoveHs kept them as graph atoms, which
    deprotonation could not remove."""

    records = list(dimorphite_dl.LoadSMIFile(StringIO(smiles + "\n")))

    assert records == [{"smiles": canonical_smiles(expected), "data": []}]


def test_deuterated_acid_is_deprotonated(
    canonical_smiles: Callable[[str], str],
) -> None:
    """Checks that an O-D acid ionizes at very basic pH. It kept the neutral
    acid with a "No valid deprotonated state" warning, because the D atom
    stayed bonded to the charged oxygen."""

    output = dimorphite_dl.run(
        smiles="[2H]C([2H])([2H])C(=O)O[2H]",
        min_ph=1e7,
        max_ph=1e7,
        return_as_list=True,
    )

    assert output is not None
    assert [line.split()[0] for line in output] == [
        canonical_smiles("[2H]C([2H])([2H])C(=O)[O-]")
    ]


def test_neutralize_mol_raises_when_a_rule_never_converges(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Checks that a rule whose product still matches its reactant raises
    with the input SMILES. The loop had no pass limit, so one such rule
    hung the whole batch."""

    real_reaction_from_smarts = dimorphite_dl.AllChem.ReactionFromSmarts

    def self_regenerating(smarts: str) -> rdChemReactions.ChemicalReaction:
        """Swaps the thiolate rule for one that leaves S- charged.

        Args:
            smarts: The reaction SMARTS neutralize_mol asked for.

        Returns:
            The compiled reaction, broken only for the thiolate rule.
        """

        if smarts.startswith("[Sv1-1:1]>>"):
            return real_reaction_from_smarts("[Sv1-1:1]>>[S-1:1]")
        return real_reaction_from_smarts(smarts)

    monkeypatch.setattr(dimorphite_dl.AllChem, "ReactionFromSmarts", self_regenerating)

    with pytest.raises(RuntimeError, match="did not converge") as info:
        dimorphite_dl.UtilFuncs.neutralize_mol(Chem.MolFromSmiles("CC[S-]"))
    assert "CC[S-]" in str(info.value)
    assert "[Sv1-1:1]>>[S-1:1]" in str(info.value)


def test_neutralize_mol_enumerates_one_product(
    monkeypatch: pytest.MonkeyPatch, canonical_smiles: Callable[[str], str]
) -> None:
    """Checks that each pass asks RDKit for a single product. Only the first
    is used, but the default enumerated up to 1000 per pass."""

    real_reaction_from_smarts = dimorphite_dl.AllChem.ReactionFromSmarts
    requested: List[int] = []

    class RecordingReaction:
        """Wraps a reaction to record how many products each call asks for,
        since RDKit's C++ methods cannot be patched directly."""

        def __init__(self, rxn: rdChemReactions.ChemicalReaction) -> None:
            """Stores the real reaction.

            Args:
                rxn: The compiled reaction to delegate to.
            """

            self.rxn = rxn

        def RunReactants(
            self, reactants: Tuple[Chem.Mol, ...], maxProducts: int = 1000
        ) -> Tuple[Tuple[Chem.Mol, ...], ...]:
            """Records maxProducts, then runs the real reaction.

            Args:
                reactants: Passed through unchanged.
                maxProducts: The product limit to record.

            Returns:
                The real reaction's products.
            """

            requested.append(maxProducts)
            return self.rxn.RunReactants(reactants, maxProducts)

    monkeypatch.setattr(
        dimorphite_dl.AllChem,
        "ReactionFromSmarts",
        lambda smarts: RecordingReaction(real_reaction_from_smarts(smarts)),
    )

    mol = dimorphite_dl.UtilFuncs.neutralize_mol(
        Chem.MolFromSmiles("[O-]C(=O)CC(=O)[O-]")
    )

    assert mol is not None
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == canonical_smiles("OC(=O)CC(=O)O")
    assert requested != [] and all(n == 1 for n in requested), requested


def test_neutralization_failure_names_its_cause(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Checks that a molecule rejected after neutralization is reported as
    such. Parse, neutralization, and RemoveHs failures all printed the same
    "poorly formed SMILES" warning."""

    monkeypatch.setattr(
        dimorphite_dl.UtilFuncs, "neutralize_mol", staticmethod(lambda mol: None)
    )

    assert list(dimorphite_dl.LoadSMIFile(StringIO("CCCN name\n"))) == []
    err = capsys.readouterr().err
    assert "sanitization failed after neutralizing charges" in err, err
    assert "CCCN name" in err, err


def test_reparse_failure_is_reported_as_unprotonated(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Checks that when Protonate cannot re-read LoadSMIFile's canonical
    SMILES, stderr says the line was written unprotonated and names it. It
    printed only "ERROR:" and the SMILES, with no hint of what happened."""

    real_convert = dimorphite_dl.UtilFuncs.convert_smiles_str_to_mol
    calls: List[str] = []

    def fail_on_reparse(smiles_str: str) -> Optional[Chem.Mol]:
        """Parses the input line, then fails Protonate's second parse.

        Args:
            smiles_str: The SMILES string to parse.

        Returns:
            The parsed Mol on the first call, then None.
        """

        calls.append(smiles_str)
        return real_convert(smiles_str) if len(calls) == 1 else None

    monkeypatch.setattr(
        dimorphite_dl.UtilFuncs,
        "convert_smiles_str_to_mol",
        staticmethod(fail_on_reparse),
    )

    output = dimorphite_dl.run(smiles="CCC(=O)O acid", return_as_list=True)

    assert output == ["CCC(=O)O\tacid"]
    err = capsys.readouterr().err
    assert "writing it unprotonated" in err, err
    assert "CCC(=O)O\tacid" in err, err


def test_caller_supplied_handle_is_left_open() -> None:
    """Checks that a file object passed in by the caller is not closed at
    EOF or by main(). The loader does not own it, and closing it broke
    callers that reuse the handle afterward."""

    handle = StringIO("CCCN\n")
    with dimorphite_dl.LoadSMIFile(handle) as loader:
        assert list(loader) == [{"smiles": "CCCN", "data": []}]
    assert not handle.closed

    handle = StringIO("CCCN\n")
    assert dimorphite_dl.run(smiles_file=handle, return_as_list=True) is not None
    assert not handle.closed


def test_defaults_agree_everywhere() -> None:
    """Checks that the command line, clean_args, and the substructure loader
    share one set of pH and precision defaults. Each had its own copy of
    the numbers, free to drift apart."""

    expected = {
        "min_ph": dimorphite_dl.DEFAULT_MIN_PH,
        "max_ph": dimorphite_dl.DEFAULT_MAX_PH,
        "pka_precision": dimorphite_dl.DEFAULT_PKA_PRECISION,
    }
    # The values the README documents.
    assert list(expected.values()) == [6.4, 8.4, 1.0]

    parser = dimorphite_dl.ArgParseFuncs.get_args()
    cleaned = dimorphite_dl.ArgParseFuncs.clean_args({"smiles": "C"})
    loader_params = inspect.signature(
        dimorphite_dl.ProtSubstructFuncs.load_protonation_substructs_calc_state_for_ph
    ).parameters

    for key, value in expected.items():
        assert parser.get_default(key) == value, key
        assert cleaned[key] == value, key
    assert loader_params["min_ph"].default == expected["min_ph"]
    assert loader_params["max_ph"].default == expected["max_ph"]
    assert loader_params["pka_std_range"].default == expected["pka_precision"]
