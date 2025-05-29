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

Dimorphite-DL adds hydrogen atoms to molecular representations, as appropriate
for a user-specified pH range. It is a fast, accurate, accessible, and modular
open-source program for enumerating small-molecule ionization states.

Users can provide SMILES strings from the command line or via an .smi file.

## Installation

You can install the latest released version on [PyPI](https://pypi.org/) using the following command.

```bash
pip install dimorphite_dl
```

Or you can install the latest development version from the `main` branch on [GitHub](https://github.com/durrantlab/dimorphite_dl) using

```bash
pip install https://github.com/durrantlab/dimorphite_dl.git
```

## Usage

```
usage: dimorphite_dl [-h] [--ph_min MIN] [--ph_max MAX]
                        [--precision PRE] [--smiles SMI]
                        [--smiles_file FILE] [--output_file FILE]
                        [--label_states] [--test]

Dimorphite 1.2.5: Creates models of appropriately protonated small moleucles.
Apache 2.0 License. Copyright 2020 Jacob D. Durrant.

optional arguments:
  -h, --help           show this help message and exit
  --ph_min MIN         minimum pH to consider (default: 6.4)
  --ph_max MAX         maximum pH to consider (default: 8.4)
  --precision PRE  pKa precision factor (number of standard devations,
                       default: 1.0)
  --smiles SMI         SMILES string to protonate
  --smiles_file FILE   file that contains SMILES strings to protonate
  --output_file FILE   output file to write protonated SMILES (optional)
  --label_states       label protonated SMILES with target state (i.e.,
                       "DEPROTONATED", "PROTONATED", or "BOTH").
  --test               run unit tests (for debugging)
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

## Examples

```bash
dimorphite_dl --smiles_file sample_molecules.smi
dimorphite_dl --smiles "CCC(=O)O" --ph_min -3.0 --ph_max -2.0
dimorphite_dl --smiles "CCCN" --ph_min -3.0 --ph_max -2.0 --output_file output.smi
dimorphite_dl --smiles_file sample_molecules.smi --precision 2.0 --label_states
dimorphite_dl --test
```

## Advanced

It is also possible to access Dimorphite-DL from another Python script, rather
than from the command line. Here's an example:

```python
from rdkit import Chem
import dimorphite_dl

# Using the dimorphite_dl.run() function, you can run Dimorphite-DL exactly as
# you would from the command line. Here's an example:
dimorphite_dl.cli.run(
   smiles="CCCN",
   ph_min=-3.0,
   ph_max=-2.0,
   output_file="output.smi"
)
print("Output of first test saved to output.smi...")

# Using the dimorphite_dl.run_with_mol_list() function, you can also pass a
# list of RDKit Mol objects. The first argument is always the list.
smiles = ["C[C@](F)(Br)CC(O)=O", "CCCCCN"]
mols = [Chem.MolFromSmiles(s) for s in smiles]
for i, mol in enumerate(mols):
    mol.SetProp("msg","Orig SMILES: " + smiles[i])

protonated_mols = dimorphite_dl.cli.run_with_mol_list(
    mols,
    ph_min=5.0,
    ph_max=9.0,
)
print([Chem.MolToSmiles(m) for m in protonated_mols])

# Note that properties are preserved.
print([m.GetProp("msg") for m in protonated_mols])
```

## Caveats

Dimorphite-DL deprotonates indoles and pyrroles around pH 14.5. But these
substructures can also be protonated around pH -3.5. Dimorphite does not
perform the protonation.

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

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi: [10.1186/s13321-019-0336-9](https://doi.org/10.1186/s13321-019-0336-9).

## License

This project is released under the Apache-2.0 License as specified in `LICENSE.md`.
