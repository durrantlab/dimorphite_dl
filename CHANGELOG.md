# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

## [2.1.0] - 2026-10-02

### Changed

- When `max_variants` truncates the output, the most probable states are kept, scored by the Henderson-Hasselbalch fraction at the middle of the pH range.
  Previously the first states generated were kept, and enumeration stopped there, so every later site stayed in its input state.
- Console logging from `enable_logging` and `--log_level` now goes to stderr, so stdout carries only SMILES.
  The `stdout_set` argument keeps its name.
- Output order no longer depends on `PYTHONHASHSEED`; variants come out in enumeration order.
- `site_substructures.smarts` now carries narrow entries for simple aliphatic amines, unactivated phenols, alkanethiols, primary sulfonamides, C-substituted amidines, N-alkyl anilines, simple carboxamides, tetrazoles, and the nucleobases, each placed ahead of the broader entry it is drawn from.
  Compounds whose pKa lies well outside the pH range now resolve to a single state rather than both: phenol, ethanethiol, benzenesulfonamide, methanesulfonamide, methylamine, benzylamine, ethanolamine, piperidine, triethylamine, quinuclidine, propranolol, 2,2,2-trifluoroethylamine, acetamidine, benzamidine, N-methylaniline, benzamide, acetanilide, peptide backbone N-H, and 5-methyltetrazole, among others.
  Amides and anilides are no longer deprotonated at physiological pH, so with `max_variants=1` a molecule bearing both an amide and a genuinely acidic group (e.g., a tetrazole) is ionized at the acidic group.
  This changes the output for existing inputs, so pipelines that depend on the previous enumeration should be checked.
  The means and standard deviations fitted from the training data are unchanged; each added entry takes its values from published measurements, cited in the substructure file.
- Measured against the external test sets of Baltruschat M, Czodrowski P (2020) F1000Research 9:113, scoring the annotated ionization site of 399 of their 403 compounds, the share returning only the dominant species at pH 7.4 rose from 48% to 63%.
  The compounds whose output omits the dominant species are the same fifteen before and after, as are the fifteen whose site is never enumerated at any pH.
- The CLI now appends each input's identifier to every variant it produces, as `SMILES,identifier`.
  Names in the input file were dropped, so variants could not be matched back to the molecule they came from.
  Scripts that parse CLI output should expect the comma-separated field; from Python, pass `label_identifiers=True` to `protonate_smiles` for the same output.
- An identifier may now contain whitespace.
  Only the second whitespace-separated field was kept, so `methyl phosphate` was truncated to `methyl`, and a list item with more than two fields was rejected rather than read.

### Fixed

- `label_states` and `--label_states` had no effect.
  Each output now lists the target state of every pKa site, in detection order.
- Neutralization of charged inputs.
    - Charge-separated sulfoxides, sulfones, phosphine oxides, and phosphates (e.g., `C[S+](C)[O-]`) are collapsed to S=O and P=O instead of protonated to hydroxysulfonium and hydroxyphosphonium cations.
    - Thiolates (`[S-]`) are protonated, so thiol and thioic acid sites are found.
    - Amidinium, guanidinium, imidazolium, and pyrazolium cations drawn with the charge on a nitrogen without H (e.g., `C[n+]1cc[nH]c1`) are neutralized.
    - N-oxides and nitrones keep their O- instead of becoming N-hydroxy cations.
    - Nitrate gets one proton, not two, so it is deprotonated again at physiological pH.
    - The terminal N- of a diazo group is left alone, so diazo compounds are no longer permanent +1 cations.
    - A rule that keeps matching its own product raises an error instead of looping forever.
- Site detection.
    - Two matches of one pattern can no longer claim the same site.
      A primary amide nitrogen was listed and charged twice, and phosphoric acid lost all three protons.
    - A carbon used by one group as its attachment point no longer hides another group's site.
      The amine of `NCP(=O)(O)O` was never protonated.
    - A nitro group no longer hides the site on the atom it is attached to.
      The N-H of N-nitrourethane (pKa 3.28) was never deprotonated, because the nitro pattern matched that nitrogen as context and protected it.
      The pattern now requires the substituent without matching it, so both sites are found on the one molecule.
- Nitroguanidines were returned as cations.
  A nitro group all but removes guanidine basicity (nitroguanidine has a pKaH of -0.98, against 13.6 for guanidine), but the general guanidine entry assigned them 12.03.
  Nitroarginine and related compounds are affected.
- N-nitroimines were returned in both states.
  The neonicotinoid insecticides carry the nitro group on the imine nitrogen, where it costs about twelve units of basicity (imidacloprid's conjugate acid is 1.56), but the general amidine entry assigned them 10.04.
  Thiamethoxam, imidacloprid, clothianidin, nitenpyram, and dinotefuran are affected.
- Nucleobases were returned in many protonation states at physiological pH.
  No entry covered their ring nitrogens, so each fell to the generic aromatic-nitrogen rules, whose windows reach the default pH range: ATP, GMP, and GTP gave 16 variants, and CMP, UMP, and dTMP gave 4, where each should give 2.
  At low pH the same rules protonated every ring nitrogen and the exocyclic amine together, returning adenine as a tetracation where only N1 is protonated.
  New entries for 9-substituted adenine and guanine and 1-substituted cytosine and uracil (including thymine) now assign one site per measured pKa: adenine N1 (pKaH 3.99), guanine N7 (pKaH 2.61) and N1-H (pKa 9.45), cytosine N3 (pKaH 4.33), and uracil N3-H (pKa 9.60).
  The remaining ring nitrogens and the exocyclic amines are claimed so that no generic rule reaches them.
  Nucleosides, nucleotides, and cofactors such as NAD, CoA, and SAM are affected, so pipelines that depend on the previous enumeration should be checked.
  The values come from published measurements and the IUPAC Digitized pKa Dataset, cited in the substructure file.
  Free bases (N9-H or N1-H), 5-halo pyrimidines, fused analogs, and hypoxanthine keep the generic rules.
- A protonation state that cannot be kekulized (e.g., a +1 charge on a bridgehead aromatic nitrogen) now keeps the site's previous state.
  Such variants were dropped as invalid SMILES, and when that was the site's only state the molecule vanished from the output.
  A molecule whose every variant is rejected now falls back to its input.
- Deuterium and tritium on O, N, or S are treated as ordinary hydrogen, so deuterated acids can be deprotonated.
  Labels on carbon are kept.
- A UTF-8 byte order mark at the start of an input file is ignored instead of corrupting the first SMILES.
- The CLI refuses an `--output_file` that is the input file, which it previously truncated before reading.
  Invalid arguments no longer truncate an existing output file, and the output file is closed properly.
- A SMILES containing `/` was read as a file path, and failed with "File not found", whenever the text before the last slash named a directory in the working directory (e.g., `CC/C=C(C)C` next to a directory named `CC`).
  A string is now read as a path only if it exists or ends in a recognized extension.

## [2.0.2] - 2025-08-11

### Added

- Turn on and control logging through the CLI.
- `colorize` keyword argument for `enable_logging` for logs to not use ANSI color codes.

### Fixed

- Determining if a provided input string was a SMILES or path to file.
  `CCC(C)=C(Cl)C/C(I)=C(\C)F` was incorrectly classified as a file.

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
