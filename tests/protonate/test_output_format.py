"""Checks the labels and ordering of protonate_smiles output."""

import os
import subprocess
import sys

import pytest

from dimorphite_dl import protonate_smiles

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

VERY_ACIDIC_PH = -10000000.0
VERY_BASIC_PH = 10000000.0
CARBOXYL_PKA = 3.456652971502591


def split_labeled_line(line: str) -> tuple[str, list[str]]:
    """Splits one labeled output line, because label_states appends the
    tab-joined states after a comma.

    Args:
        line: An output line from protonate_smiles with label_states=True.

    Returns:
        The SMILES and the list of state labels.
    """

    smiles, _, states = line.partition(",")
    return smiles, states.split("\t") if states else []


@pytest.mark.parametrize(
    ("ph", "expected_state"),
    [
        (VERY_ACIDIC_PH, "PROTONATED"),
        (VERY_BASIC_PH, "DEPROTONATED"),
        (CARBOXYL_PKA, "BOTH"),
    ],
)
def test_label_states_reports_site_state(ph: float, expected_state: str) -> None:
    """Checks that label_states adds the target state. It read a
    nonexistent site attribute, so every label was empty."""

    output = protonate_smiles(
        "CCC(=O)O", ph_min=ph, ph_max=ph, precision=0.5, label_states=True
    )

    assert len(output) > 0
    for line in output:
        _, states = split_labeled_line(line)
        assert states == [expected_state], output


def test_label_states_lists_every_pka() -> None:
    """Checks that a site with two pKa values (phosphate) and a second site
    each contribute their own label."""

    output = protonate_smiles(
        "NCCOP(=O)(O)O",
        ph_min=VERY_BASIC_PH,
        ph_max=VERY_BASIC_PH,
        precision=0.5,
        label_states=True,
    )

    assert len(output) == 1
    _, states = split_labeled_line(output[0])
    assert states == ["DEPROTONATED"] * 3, output


def test_no_labels_without_label_states() -> None:
    """Checks that the default output stays a bare SMILES."""

    output = protonate_smiles(
        "CCC(=O)O", ph_min=VERY_BASIC_PH, ph_max=VERY_BASIC_PH, precision=0.5
    )

    assert len(output) == 1
    assert "," not in output[0] and "\t" not in output[0]


def test_output_order_is_independent_of_hash_seed() -> None:
    """Checks that variant order does not change between runs. Variants were
    deduplicated through a set, so their order followed string hashing,
    which Python randomizes per process."""

    script = (
        "from dimorphite_dl import protonate_smiles\n"
        "for line in protonate_smiles('OC(=O)CCC(N)C(=O)O', ph_min=3.5, "
        "ph_max=3.5, precision=1.0, label_states=True):\n"
        "    print(line)\n"
    )

    outputs = []
    for seed in ["0", "1", "2", "3", "4"]:
        env = dict(os.environ)
        env["PYTHONHASHSEED"] = seed
        # Logging would add timestamped lines that differ between runs.
        env.pop("DIMORPHITE_DL_LOG", None)
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=PROJECT_ROOT,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=True,
        )
        outputs.append(result.stdout)

    assert len(outputs[0].splitlines()) > 1, outputs[0]
    assert all(output == outputs[0] for output in outputs), outputs


def test_states_stay_in_their_column_without_identifier() -> None:
    """Checks that the states column does not move when an input line has
    no name. The empty identifier was dropped, so the states landed in the
    identifier field."""

    output = protonate_smiles(
        ["CCC(=O)O", "CCC(=O)O acid"],
        ph_min=VERY_BASIC_PH,
        ph_max=VERY_BASIC_PH,
        precision=0.5,
        label_identifiers=True,
        label_states=True,
    )

    assert [line.split(",")[1:] for line in output] == [
        ["", "DEPROTONATED"],
        ["acid", "DEPROTONATED"],
    ], output
