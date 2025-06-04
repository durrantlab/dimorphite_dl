<h1 align="center">dimorphite_dl</h1>

<h4 align="center">Adds hydrogen atoms to molecular representations as specified by pH</h4>

<p align="center">
    <a href="https://github.com/durrantlab/dimorphite_dl/actions/workflows/tests.yml">
        <img src="https://github.com/durrantlab/dimorphite_dl/actions/workflows/tests.yml/badge.svg" alt="Build Status ">
    </a>
    <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/dimorphite_dl">
    <a href="https://codecov.io/gh/durrantlab/dimorphite_dl">
        <img src="https://codecov.io/gh/durrantlab/dimorphite_dl/branch/main/graph/badge.svg" alt="codecov">
    </a>
    <a href="https://github.com/durrantlab/dimorphite_dl/releases">
        <img src="https://img.shields.io/github/v/release/durrantlab/dimorphite_dl" alt="GitHub release (latest by date)">
    </a>
    <a href="https://pypistats.org/packages/dimorphite-dl">
        <img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/dimorphite-dl?style=flat&label=PyPI%20downloads&link=https%3A%2F%2Fpypistats.org%2Fpackages%2Fdimorphite-dl">
    </a>
    <a href="https://github.com/durrantlab/dimorphite_dl/blob/main/LICENSE" target="_blank">
        <img src="https://img.shields.io/github/license/durrantlab/dimorphite_dl" alt="License">
    </a>
    <a href="https://github.com/durrantlab/dimorphite_dl/" target="_blank">
        <img src="https://img.shields.io/github/repo-size/durrantlab/dimorphite_dl" alt="GitHub repo size">
    </a>
    <a href="https://doi.org/10.5281/zenodo.15486131">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.15486131.svg" alt="DOI">
    </a>
    <a href="https://archive.softwareheritage.org/browse/origin/?origin_url=https://doi.org/10.5281/zenodo.15486131">
        <img src="https://archive.softwareheritage.org/badge/origin/https://doi.org/10.5281/zenodo.15486131/" alt="Archived | https://doi.org/10.5281/zenodo.15486131"/>
    </a>
</p>

Dimorphite-DL is a fast, accurate, accessible, and modular open-source program designed for enumerating small-molecule ionization states.
It specifically adds or removes hydrogen atoms from molecular representations to achieve the appropriate protonation state for a user-specified pH range.

Accurate protonation states are crucial in cheminformatics and computational drug discovery, as a molecule's ionization state significantly impacts its physicochemical properties, biological activity, and interactions with targets.
Dimorphite-DL addresses this by providing a robust solution for preparing molecules for various downstream applications like docking, molecular dynamics, and virtual screening.

## Installation

You can install the latest released version on [PyPI](https://pypi.org/project/dimorphite-dl/) using the following command.

```bash
pip install dimorphite_dl
```

Or you can install the latest development version from the `main` branch on [GitHub](https://github.com/durrantlab/dimorphite_dl) using

```bash
pip install https://github.com/durrantlab/dimorphite_dl.git
```

## Usage

### CLI

The command-line interface provides straightforward access to Dimorphite-DL's functionalities.

```bash
dimorphite_dl [-h] [--ph_min MIN] [--ph_max MAX] [--precision PRE] [--output_file FILE] [--max_variants MXV] [--label_states] [--silent] SMI

dimorphite_dl v1.2.5
```

**Positional Arguments:**

- `SMI`: SMILES string or path to a file containing SMILES strings to protonate.

**Options:**

- `--ph_min MIN`: Minimum pH to consider (default: 6.4).
- `--ph_max MAX`: Maximum pH to consider (default: 8.4).
- `--precision PRE`: pKa precision factor, representing the number of standard deviations from the mean pKa to consider when determining ionization states (default: 1.0).
- `--output_file FILE`: Optional path to a file to write the protonated SMILES results.
- `--max_variants MXV`: Limits the number of protonation variants generated per input compound (default: 128).
- `--label_states`: If set, output SMILES will be labeled with their target ionization state ("DEPROTONATED", "PROTONATED", or "BOTH").
- `--silent`: Suppresses all messages printed to the screen.

#### Examples

**Protonate molecules from a file:**

```bash
dimorphite_dl sample_molecules.smi
```

**Protonate a single SMILES string within a specific pH range:**

```bash
dimorphite_dl --ph_min -3.0 --ph_max -2.0 "CCC(=O)O"
```

**Protonate a SMILES string and save output to a file:**

```bash
dimorphite_dl --ph_min -3.0 --ph_max -2.0 --output_file output.smi "CCCN"
```

**Protonate molecules from a file with increased pKa precision and state labels:**

```bash
dimorphite_dl --precision 2.0 --label_states sample_molecules.smi
```

### Scripting

Dimorphite-DL can be easily integrated into your Python scripts.
The primary function for this is `protonate_smiles` from `dimorphite_dl.protonate`.

```python
from dimorphite_dl import protonate_smiles

# Protonate a single SMILES string with custom pH range and precision
protonated_mol_1: list[str] = protonate_smiles(
    "CCC(=O)O", ph_min=6.8, ph_max=7.9, precision=0.5
)
print(f"Protonated 'CCC(=O)O': {protonated_mol_1}")

# Protonate a list of SMILES strings
protonated_mol_list: list[str] = protonate_smiles(["CCC(=O)O", "CCCN"])
print(f"Protonated list: {protonated_mol_list}")

# Protonate molecules from a SMILES file
# Make sure '~/example.smi' exists and contains SMILES strings
# protonated_from_file: list[str] = protonate_smiles("~/example.smi")
# print(f"Protonated from file: {protonated_from_file}")

# Example with labeling states and limiting variants
protonated_labeled: list[str] = protonate_smiles(
    "C1CCCCC1C(=O)O", ph_min=7.0, ph_max=7.4, label_states=True, max_variants=5
)
print(f"Protonated with labels: {protonated_labeled}")
```

## Known issues

Dimorphite_dl is designed to handle the vast majority of ionizable functional groups accurately, but there are some edge cases where the current SMARTS patterns and pKa assignments may not behave as expected.
The following are known limitations that users should be aware of when working with specific molecular substructures:

- **Tertiary Amides**: Tertiary amides (e.g., N-acetylpiperidine `CC(=O)N1CCCCC1`) are incorrectly treated as basic amines (pKa ~8) instead of neutral species because current amide SMARTS patterns require an N-H bond.
- **Indoles and Pyrroles**: These heterocycles are correctly deprotonated around pH 14.5 but are not protonated at very low pH (~-3.5) where they would be expected to protonate under extremely acidic conditions.

## Development

We use [pixi](https://pixi.sh/latest/) to manage Python environments and simplify the developer workflow.
Once you have [pixi](https://pixi.sh/latest/) installed, move into `dimorphite_dl` directory (e.g., `cd dimorphite_dl`) and install the environment using the command

```bash
pixi install
```

Now you can activate the new virtual environment using

```sh
pixi shell
```

## Citation

If you use Dimorphite-DL in your research, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An open-source program for enumerating the ionization states of drug-like small
molecules. *J Cheminform 11*:14. doi: [10.1186/s13321-019-0336-9](https://doi.org/10.1186/s13321-019-0336-9).

## License

This project is released under the Apache-2.0 License as specified in `LICENSE.md`.
