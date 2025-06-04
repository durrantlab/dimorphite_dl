# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `colorize` keyword argument for `enable_logging` for logs to not use ANSI color codes.

## [2.0.1] - 2025-06-03

### Changed

- Rearranged `__init__.py` imports to mainly have `from dimorphite_dl import protonate_smiles`.

### Fixed

- Circular import of SMARTS.

## [2.0.0] - 2025-06-01

### Changed

- Fallback mechanism now uses the previous successful site protonation.
  In previous versions, sometimes only the last successful protonation site type was returned.
  If the third phosphate protonation failed, then it would fallback to the last successful protonation before the first phosphate.
  Now, we would return the second phosphate protonation.
- Major refactor of practically everything.

## [1.2.5] - 2025-05-21

### Changed

- Major reorganization of the original `dimorphite_dl.py` file into Python modules under the package name `dimorphite_dl`. No code logic has been change, just refactored.

## [1.2.4]

### Added

- Added test cases for ATP and NAD.

### Changed

- Dimorphite-DL now better protonates compounds with polyphosphate chains
  (e.g., ATP). See `site_substructures.smarts` for the rationale behind the
  added pKa values.
- `site_substructures.smarts` now allows comments (lines that start with `#`).
- Improved suport for the `--silent` option.
- Reformatted code per the [*Black* Python code formatter](https://github.com/psf/black).

### Fixed

- Fixed a bug that affected how Dimorphite-DL deals with new protonation
    states that yield invalid SMILES strings.
    - Previously, it simply returned the original input SMILES in these rare
    cases (better than nothing). Now, it instead returns the last valid SMILES
    produced, not necessarily the original SMILES.
    - Consider `O=C(O)N1C=CC=C1` at pH 3.5 as an example.
        - Dimorphite-DL first deprotonates the carboxyl group, producing
      `O=C([O-])n1cccc1` (a valid SMILES).
        - It then attempts to protonate the aromatic nitrogen, producing
      `O=C([O-])[n+]1cccc1`, an invalid SMILES.
        - Previously, it would output the original SMILES, `O=C(O)N1C=CC=C1`. Now
      it outputs the last valid SMILES, `O=C([O-])n1cccc1`.

## [1.2.3]

### Added

- Added "silent" option to suppress all output.
- Added code to suppress unnecessary RDKit warnings.

### Changed

- Updated protonation of nitrogen, oxygen, and sulfur atoms to be compatible
  with the latest version of RDKit, which broke backwards compatibility.
- Updated copyright to 2020.

## [1.2.2]

### Added

- Added a new parameter to limit the number of variants per compound
  (`--max_variants`). The default is 128.

## [1.2.1]

### Fixed

- Corrected a bug that rarely misprotonated/deprotonated compounds with
  multiple ionization sites (e.g., producing a carbanion).

## [1.2.0]

### Fixed

- Corrected a bug that led Dimorphite-DL to sometimes produce output molecules
  that are non-physical.
- Corrected a bug that gave incorrect protonation states for rare molecules
  (aromatic rings with nitrogens that are protonated when electrically
  neutral, e.g. pyridin-4(1H)-one).
- `run_with_mol_list()` now preserves non-string properties.
- `run_with_mol_list()` throws a warning if it cannot process a molecule,
  rather than terminating the program with an error.

## [1.1.0]

### Added

- Dimorphite-DL now distinguishes between indoles/pyrroles and
  Aromatic_nitrogen_protonated.
- It is now possible to call Dimorphite-DL from another Python script, in
  addition to the command line. See the `README.md` file for instructions.

## [1.0.0]

The original version described in:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.
