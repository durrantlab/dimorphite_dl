Dimorphite-DL 1.1
=================

What is it?
-----------

Dimorphite-DL adds hydrogen atoms to molecular representations, as appropriate
for a user-specified pH range. It is a fast, accurate, accessible, and modular
open-source program for enumerating small-molecule ionization states.

Users can provide SMILES strings from the command line or via an .smi file.

Citation
--------

If you use Dimorphite-DL in your research, please cite:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.

Licensing
---------

Protonation is released under the Apache 2.0 license. See LICENSE.txt for
details.

Usage
-----

```
usage: dimorphite_dl.py [-h] [--min_ph MIN] [--max_ph MAX]
                        [--pka_precision PRE] [--max_variants MXV]
                        [--smiles SMI] [--smiles_file FILE]
                        [--output_file FILE] [--label_states]

Dimorphite 1.1: Creates models of appropriately protonated small moleucles.
Apache 2.0 License. Copyright 2018 Jacob D. Durrant.

optional arguments:
  -h, --help           show this help message and exit
  --min_ph MIN         minimum pH to consider (default: 6.4)
  --max_ph MAX         maximum pH to consider (default: 8.4)
  --pka_precision PRE  pKa precision factor (number of standard devations,
                       default: 1.0)
  --max_variants MXV   limit number of variants per input compound, keeping
                       the most probable (default: 128)
  --smiles SMI         SMILES string to protonate
  --smiles_file FILE   file that contains SMILES strings to protonate
  --output_file FILE   output file to write protonated SMILES (optional)
  --label_states       label protonated SMILES with target state (i.e.,
                       "DEPROTONATED", "PROTONATED", or "BOTH").
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

Examples
--------

```
  python dimorphite_dl.py --smiles_file sample_molecules.smi
  python dimorphite_dl.py --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
  python dimorphite_dl.py --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
  python dimorphite_dl.py --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states
```

Advanced Usage
--------------

It is also possible to access Dimorphite-DL from another Python script, rather
than from the command line. Here's an example:

```python
from rdkit import Chem
import dimorphite_dl

# Using the dimorphite_dl.run() function, you can run Dimorphite-DL exactly as
# you would from the command line. Here's an example:
dimorphite_dl.run(
   smiles="CCCN",
   min_ph=-3.0,
   max_ph=-2.0,
   output_file="output.smi"
)
print("Output of first test saved to output.smi...")

# Using the dimorphite_dl.run_with_mol_list() function, you can also pass a
# list of RDKit Mol objects. The first argument is always the list.
mols = [Chem.MolFromSmiles(s) for s in ["C[C@](F)(Br)CC(O)=O", "CCCCCN"]]
protonated_mols = dimorphite_dl.run_with_mol_list(
    mols,
    min_ph=5.0,
    max_ph=9.0,
)
print([Chem.MolToSmiles(m) for m in protonated_mols])

# Each input can yield several protonated mols (or none, if it fails), so each
# output records the position of the input it came from.
print([m.GetIntProp("dimorphite_input_index") for m in protonated_mols])
```

Testing
-------

The tests use pytest. From the project root, run:

```
pytest tests/
```

Caveats
-------

Dimorphite-DL deprotonates indoles and pyrroles around pH 14.5. But these
substructures can also be protonated around pH -3.5. Dimorphite does not
perform the protonation.

The following limitations are part of the published substructure library
(`site_substructures.smarts`). They are documented here rather than fixed
so that results remain consistent with the original publication.

* Tertiary amides, sulfonamides, carbamates, and ureas are protonated as if
  they were aliphatic amines (pKa ~8.2). The amide and sulfonamide
  substructures require an N-H, so a nitrogen without one falls through to the
  general amine substructure, which matches any trivalent nitrogen bonded to
  an aliphatic carbon. For example, `CN(C)C=O` becomes `C[NH+](C)C=O` at pH
  2, and `CC(=O)N1CCN(C)CC1` becomes doubly protonated. At the default pH
  range the correct neutral form is still produced, but the protonated variant
  is ranked as more probable, so it is the one kept when `--max_variants` is
  small.
* Acidic aromatic N-H groups other than indoles and pyrroles are never
  deprotonated. Tetrazoles (pKa ~4.9), for example, remain neutral at every
  pH, and cationic tetrazolium variants may be produced instead. Aromatic
  imides such as uracil and 5-fluorouracil are also not recognized, because
  the ringed-imide substructures match only non-aromatic rings. Molecules with
  aromatic N-H groups (imidazoles, pyrazoles, triazoles, pyridones, etc.) may
  also trigger "No valid protonated state" warnings; these do not affect the
  output.
* N-alkyl anilines whose alkyl carbon bears exactly one hydrogen (e.g.,
  N-isopropyl and N-cyclohexyl anilines) are treated as aliphatic amines
  rather than anilines. In SMARTS, `[!H]` means "not exactly one attached
  hydrogen," not "not a hydrogen atom," so the secondary and tertiary aniline
  substructures do not match these compounds. At pH 7.4, for example,
  `CC(C)Nc1ccccc1` produces a protonated variant, but `CCNc1ccccc1` does not.
* Phosphorus atoms in di- and triphosphate linkages that carry only one OH
  group are never deprotonated. This includes the alpha and beta phosphates
  of ATP and both phosphates of NAD. The phosphate diester substructure
  requires both ester oxygens to be bonded to carbon, nitrogen, or a halogen,
  so a P-O-P bridge is not recognized. ATP is therefore assigned a net charge
  of -2 rather than -4, even at pH 12. Charges for nucleotides, cofactors, and
  other polyphosphates should be assigned manually.

Authors and Contacts
--------------------

See the `CONTRIBUTORS.md` file for a full list of contributors. Please contact
Jacob Durrant (durrantj@pitt.edu) with any questions.
