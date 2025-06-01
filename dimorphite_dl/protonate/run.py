"""
Robust protonation processing module.

This module provides the main Protonate class for processing SMILES strings
and generating protonated variants. Each function is focused on a specific
aspect of the protonation workflow with comprehensive error handling.
"""

from typing import Generator, Iterable, Iterator

from dataclasses import dataclass

from loguru import logger
from rdkit import Chem

from dimorphite_dl.io import SMILESProcessor, SMILESRecord
from dimorphite_dl.mol import MoleculeRecord
from dimorphite_dl.protonate.change import protonate_site
from dimorphite_dl.protonate.data import PKaData
from dimorphite_dl.protonate.detect import ProtonationSiteDetector
from dimorphite_dl.protonate.site import ProtonationSite


class ProtonationError(Exception):
    """Raised when protonation processing encounters an error."""

    pass


@dataclass
class ProtonationResult:
    """Data structure for protonation results with explicit fields."""

    smiles: str
    identifier: str
    states: str = ""

    def __post_init__(self):
        """Validate result data after initialization."""
        assert isinstance(self.smiles, str)
        assert isinstance(self.identifier, str)
        assert isinstance(self.states, str)
        assert len(self.smiles) > 0  # SMILES cannot be empty

    def to_string(
        self,
        include_identifier: bool = False,
        include_states: bool = False,
        separator: str = ",",
    ) -> str:
        """
        Convert result to output string format.

        Args:
            include_identifier: Whether to include the identifier.
            include_states: Whether to include state information.
            separator: What to separate additional information with.

        Returns:
            Formatted string representation
        """
        assert isinstance(include_identifier, bool)
        assert isinstance(include_states, bool)
        assert isinstance(separator, str)

        output = self.smiles
        if include_identifier and self.identifier != "":
            output += separator + self.identifier
        if include_states and len(self.states) > 0:
            output += separator + self.states
        return output


@dataclass
class ProtonationStats:
    """Statistics tracking for protonation processing with explicit counters."""

    molecules_processed: int = 0
    total_variants_generated: int = 0
    variants_validated: int = 0
    variants_rejected: int = 0
    molecules_with_sites: int = 0
    molecules_without_sites: int = 0
    fallback_used: int = 0

    def __post_init__(self):
        """Validate statistics after initialization."""
        assert self.molecules_processed >= 0
        assert self.total_variants_generated >= 0
        assert self.variants_validated >= 0
        assert self.variants_rejected >= 0
        assert self.molecules_with_sites >= 0
        assert self.molecules_without_sites >= 0
        assert self.fallback_used >= 0


class Protonate:
    """
    Generator class for protonating SMILES strings with comprehensive error handling.

    Processes molecules one at a time, generating protonated variants based on
    pH conditions and pKa data. Includes validation, fallback mechanisms,
    and detailed statistics tracking.

    Design goals: Safety through validation, performance through batching,
    developer experience through clear error messages and statistics.
    """

    def __init__(
        self,
        smiles_input: str | Iterable[str] | Iterator[str] | None = None,
        ph_min: float = 6.4,
        ph_max: float = 8.4,
        precision: float = 1.0,
        label_identifiers: bool = False,
        label_states: bool = False,
        max_variants: int = 128,
        validate_output: bool = True,
        **smiles_processor_kwargs,
    ):
        """
        Initialize the protonation generator with explicit parameters.

        Args:
            smiles_input: SMILES string, file path, or iterable of SMILES
            ph_min: Minimum pH to consider.
            ph_max: Maximum pH to consider.
            precision: pKa precision factor.
            label_identifiers: When returning SMILES, format the string to
                include any identifier.
            label_states: When returning SMILES, format the string to include
                states.
            max_variants: Maximum number of variants per input compound (bounded)
            validate_output: Whether to validate generated SMILES (explicit)
            **smiles_processor_kwargs: Additional arguments for SMILESProcessor
        """
        # Validate all input parameters with clear bounds
        assert ph_min <= ph_max, (
            f"ph_min ({ph_min}) must be less than or equal to ph_max ({ph_max})"
        )
        assert precision >= 0.0 and precision <= 10.0, (
            f"precision must be 0-10, got: {precision}"
        )
        assert max_variants > 0 and max_variants <= 10000, (
            f"max_variants must be 1-10000, got: {max_variants}"
        )
        assert isinstance(label_identifiers, bool)
        assert isinstance(label_states, bool)
        assert isinstance(validate_output, bool)

        self.smiles_input = smiles_input
        self.ph_min = ph_min
        self.ph_max = ph_max
        self.precision = precision
        self.label_identifiers = label_identifiers
        self.label_states = label_states
        self.max_variants = max_variants
        self.validate_output = validate_output

        logger.info(
            "Initializing protonation with pH {:.1f} to {:.1f}, precision: {:.1f}, max_variants: {}",
            ph_min,
            ph_max,
            precision,
            max_variants,
        )

        self.smiles_processor = SMILESProcessor(**smiles_processor_kwargs)
        self.stats = ProtonationStats()
        self.pka_data = PKaData()
        self.site_detector = ProtonationSiteDetector(
            validate_sites=True, max_sites_per_molecule=50
        )

        # Current processing state
        self.current_results_queue: list[ProtonationResult] = []
        self._smiles_stream: Iterator[SMILESRecord] | None = None

        self._initialize_smiles_stream()

    def _initialize_smiles_stream(self) -> None:
        """
        Initialize the SMILES input stream with error handling.

        Raises:
            ProtonationError: If stream initialization fails
        """
        if self.smiles_input is None:
            raise ProtonationError("No SMILES input provided")

        try:
            self._smiles_stream = iter(self.smiles_processor.stream(self.smiles_input))
            logger.debug("Initialized SMILES stream successfully")

        except Exception as error:
            logger.error("Failed to initialize SMILES stream: {}", str(error))
            raise ProtonationError(
                f"Failed to initialize SMILES stream: {error}"
            ) from error

    def __iter__(self):
        """Return this generator object for iteration."""
        return self

    def __next__(self) -> str:
        """
        Generate the next protonated SMILES string.

        Returns:
            String containing protonated SMILES with metadata

        Raises:
            StopIteration: When no more SMILES are available to process
            ProtonationError: When processing encounters a critical error
        """
        # Return queued results first
        if len(self.current_results_queue) > 0:
            result = self.current_results_queue.pop(0)
            return result.to_string(self.label_identifiers, self.label_states)

        # Process next input SMILES if queue is empty
        next_smiles_record = self._get_next_smiles_record()
        self._process_single_smiles_record(next_smiles_record)

        # Return first result from newly populated queue
        if len(self.current_results_queue) > 0:
            result = self.current_results_queue.pop(0)
            return result.to_string(self.label_identifiers, self.label_states)
        else:
            # If no results generated, try next record
            return self.__next__()

    def _get_next_smiles_record(self) -> SMILESRecord:
        """
        Get the next SMILES record from input stream.

        Returns:
            Next SMILESRecord from input stream

        Raises:
            StopIteration: When no more input SMILES available
        """
        if self._smiles_stream is None:
            raise StopIteration("No SMILES stream available")

        try:
            smiles_record = next(self._smiles_stream)
            assert isinstance(smiles_record, SMILESRecord)
            assert len(smiles_record.smiles) > 0
            return smiles_record

        except StopIteration:
            raise StopIteration("No more SMILES to process")

    def _process_single_smiles_record(self, smiles_record: SMILESRecord) -> None:
        """
        Process one SMILES record and populate results queue.

        Args:
            smiles_record: SMILESRecord containing SMILES and metadata
        """
        assert isinstance(smiles_record, SMILESRecord)
        assert len(smiles_record.smiles) > 0

        mol_record = MoleculeRecord(smiles_record.smiles, smiles_record.identifier)
        self.stats.molecules_processed += 1

        try:
            # Find protonation sites
            mol_record, protonation_sites = self._detect_protonation_sites(mol_record)

            if mol_record.mol is None or len(protonation_sites) == 0:
                self._handle_molecule_without_sites(mol_record)
                return

            # Generate protonated variants
            protonated_molecules = self._generate_protonated_variants(
                mol_record, protonation_sites
            )

            # Convert to SMILES and validate
            smiles_strings = self._convert_molecules_to_smiles(
                list(protonated_molecules), mol_record.smiles_original
            )
            # Handle empty results with fallback
            if len(smiles_strings) == 0:
                self._add_fallback_result(
                    mol_record.smiles_original, mol_record.identifier
                )
                return

            if self.validate_output:
                smiles_strings = self._validate_generated_smiles(
                    smiles_strings, mol_record.smiles_original
                )
            # Create results with state information
            self._create_results_from_smiles(
                smiles_strings, mol_record.identifier, protonation_sites
            )

        except Exception as error:
            logger.error(
                "Error processing SMILES '{}': {}",
                mol_record.smiles_original,
                str(error),
            )
            self._add_fallback_result(mol_record.smiles_original, mol_record.identifier)

    def _detect_protonation_sites(
        self, mol_record: MoleculeRecord
    ) -> tuple[MoleculeRecord, list]:
        """
        Detect protonation sites in molecule.

        Args:
            mol_record: MoleculeRecord to analyze

        Returns:
            Tuple of (updated_mol_record, protonation_sites_list)
        """
        assert isinstance(mol_record, MoleculeRecord)

        try:
            mol_record, sites = self.site_detector.find_sites(mol_record)
            site_count = len(sites)

            if site_count > 0:
                self.stats.molecules_with_sites += 1
            else:
                logger.debug("No protonation sites found for '{}'", mol_record.smiles)

            return mol_record, sites

        except Exception as error:
            logger.warning(
                "Error detecting sites for '{}': {}", mol_record.smiles, str(error)
            )
            return mol_record, []

    def _handle_molecule_without_sites(self, mol_record: MoleculeRecord) -> None:
        """
        Handle molecules that have no protonation sites.

        Args:
            mol_record: MoleculeRecord without protonation sites
        """
        assert isinstance(mol_record, MoleculeRecord)

        logger.warning(
            "No protonation sites found for '{}'", mol_record.smiles_original
        )
        self.stats.molecules_without_sites += 1

        # Try to generate a clean molecule without explicit hydrogens
        clean_smiles = self._generate_clean_smiles_without_hydrogens(mol_record)
        if clean_smiles is not None:
            self._add_result_to_queue(clean_smiles, mol_record.identifier, "")
        else:
            self._add_fallback_result(mol_record.smiles_original, mol_record.identifier)

    def _generate_clean_smiles_without_hydrogens(
        self, mol_record: MoleculeRecord
    ) -> str | None:
        """
        Generate clean SMILES without explicit hydrogens.

        Args:
            mol_record: MoleculeRecord to process

        Returns:
            Clean SMILES string or None if generation failed
        """
        assert isinstance(mol_record, MoleculeRecord)

        try:
            if mol_record.mol is not None:
                mol_without_hs = MoleculeRecord.remove_hydrogens(mol_record.mol)
                if mol_without_hs is not None:
                    clean_smiles = Chem.MolToSmiles(
                        mol_without_hs, isomericSmiles=True, canonical=True
                    )
                    if clean_smiles and len(clean_smiles) > 0:
                        return clean_smiles
        except Exception as error:
            logger.warning(
                "Error generating clean SMILES for '{}': {}",
                mol_record.smiles,
                str(error),
            )

        return None

    def _generate_protonated_variants(
        self, mol_record: MoleculeRecord, sites: list[ProtonationSite]
    ) -> list[Chem.Mol]:
        """
        Generate protonated variants, rolling back an entire “site.name” group if any one site in that group fails.

        Args:
            mol_record: MoleculeRecord with detected sites (mol_record.mol != None)
            sites: List of ProtonationSite instances, in detection order

        Returns:
            A list of RDKit Mol objects. If an entire group (same site.name) fails, we fall back
            to the molecules that existed before that group began.
        """
        assert isinstance(mol_record, MoleculeRecord)
        assert isinstance(sites, list) and len(sites) > 0
        assert mol_record.mol is not None

        # Start with the original molecule
        mols_validated: list[Chem.Mol] = [mol_record.mol]
        max_variants = self.max_variants

        for i, site in enumerate(sites):
            # Try protonating *all* currently validated molecules at this one site
            # (no need to wrap in list(...); mols_validated is already a list)
            mols_tmp = self._protonate_single_site(
                mols_validated, site, mol_record.smiles_original, i
            )

            # Something failed and we want to exit and return previously validated
            # molecules.
            if len(mols_tmp) == 0:
                break

            # Site succeeded → keep its output
            mols_validated = mols_tmp

            # Cap at max_variants once we exceed it
            if len(mols_validated) > max_variants:
                mols_validated = mols_validated[:max_variants]
                break

        return mols_validated

    def _protonate_single_site(
        self, molecules: list[Chem.Mol], site, original_smiles: str, site_index: int
    ) -> list[Chem.Mol]:
        """
        Apply protonation to a single site across all molecules.

        Args:
            molecules: List of RDKit mol objects to protonate
            site: Protonation site to process
            original_smiles: Original SMILES for logging context
            site_index: Site index for logging context

        Returns:
            List of protonated mol objects
        """
        assert isinstance(molecules, list)
        assert len(molecules) > 0
        assert isinstance(original_smiles, str)
        assert isinstance(site_index, int)
        assert site_index >= 0

        try:
            new_molecules = protonate_site(
                molecules, site, self.ph_min, self.ph_max, self.precision
            )
            return new_molecules
        except Exception as error:
            logger.warning(
                "Failed to protonate site {} for '{}': {}",
                site_index,
                original_smiles,
                str(error),
            )
            return []

    def _convert_molecules_to_smiles(
        self, molecules: list[Chem.Mol], original_smiles: str
    ) -> list[str]:
        """
        Convert RDKit mol objects to canonical SMILES strings.

        Args:
            molecules: List of RDKit mol objects
            original_smiles: Original SMILES for logging context

        Returns:
            List of unique canonical SMILES strings
        """
        assert isinstance(molecules, list)
        assert isinstance(original_smiles, str)

        if len(molecules) == 0:
            return []

        unique_smiles = set()
        valid_molecule_count = 0

        for mol in molecules:
            if mol is None:
                continue

            valid_molecule_count += 1
            canonical_smiles = self._generate_canonical_smiles_from_mol(mol)
            if canonical_smiles is not None:
                unique_smiles.add(canonical_smiles)

        smiles_list = list(unique_smiles)
        unique_count = len(smiles_list)
        self.stats.total_variants_generated += unique_count

        logger.debug(
            "Generated {} unique SMILES from {} valid molecules for '{}'",
            unique_count,
            valid_molecule_count,
            original_smiles,
        )

        return smiles_list

    def _generate_canonical_smiles_from_mol(self, mol: Chem.Mol) -> str | None:
        """
        Generate canonical SMILES from a single mol object.

        Args:
            mol: RDKit mol object

        Returns:
            Canonical SMILES string or None if generation failed
        """
        assert mol is not None

        try:
            canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if canonical and len(canonical) > 0:
                return canonical
        except Exception as error:
            logger.debug("Error generating canonical SMILES: {}", str(error))

        return None

    def _validate_generated_smiles(
        self, smiles_list: list[str], original_smiles: str
    ) -> list[str]:
        """
        Validate list of generated SMILES strings.

        Args:
            smiles_list: List of SMILES to validate
            original_smiles: Original SMILES for logging context

        Returns:
            List of validated SMILES strings
        """
        assert isinstance(smiles_list, list)
        assert isinstance(original_smiles, str)

        if len(smiles_list) == 0:
            return []

        validated_smiles = []

        for smiles in smiles_list:
            assert isinstance(smiles, str)
            assert len(smiles) > 0

            if self._is_smiles_valid(smiles):
                validated_smiles.append(smiles)
                self.stats.variants_validated += 1
            else:
                self.stats.variants_rejected += 1
                logger.debug("Rejected invalid SMILES: '{}'", smiles)

        validated_count = len(validated_smiles)
        original_count = len(smiles_list)
        logger.debug(
            "{}/{} SMILES were valid for '{}'",
            validated_count,
            original_count,
            original_smiles,
        )

        return validated_smiles

    def _is_smiles_valid(self, smiles: str) -> bool:
        """
        Check if a SMILES string represents a valid molecule.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if SMILES is valid, False otherwise
        """
        assert isinstance(smiles, str)
        assert len(smiles) > 0

        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            return False

    def _create_results_from_smiles(
        self, smiles_list: list[str], identifier: str, sites: list
    ) -> None:
        """
        Create ProtonationResult objects from SMILES list.

        Args:
            smiles_list: List of SMILES strings
            identifier: Molecule identifier
            sites: List of protonation sites for state labeling
        """
        assert isinstance(smiles_list, list)
        assert isinstance(identifier, str)
        assert isinstance(sites, list)

        if len(smiles_list) == 0:
            return

        # Generate states string if labeling is requested
        states_string = ""
        if self.label_states and len(sites) > 0:
            states_string = self._generate_states_string_from_sites(sites)

        # Create results for each SMILES
        for smiles in smiles_list:
            assert isinstance(smiles, str)
            assert len(smiles) > 0

            result = ProtonationResult(
                smiles=smiles, identifier=identifier, states=states_string
            )
            self.current_results_queue.append(result)

    def _generate_states_string_from_sites(self, sites: list) -> str:
        """
        Generate states string from protonation sites.

        Args:
            sites: List of protonation sites

        Returns:
            Tab-separated string of site states
        """
        assert isinstance(sites, list)
        assert len(sites) > 0

        try:
            state_strings = []
            for site in sites:
                if hasattr(site, "target_state") and site.target_state:
                    state_strings.append(str(site.target_state))

            return "\t".join(state_strings)
        except Exception as error:
            logger.debug("Error generating states string: {}", str(error))
            return ""

    def _add_result_to_queue(self, smiles: str, identifier: str, states: str) -> None:
        """
        Add a single result to the processing queue.

        Args:
            smiles: SMILES string
            identifier: Molecule identifier
            states: States string
        """
        assert isinstance(smiles, str)
        assert isinstance(identifier, str)
        assert isinstance(states, str)
        assert len(smiles) > 0

        result = ProtonationResult(smiles=smiles, identifier=identifier, states=states)
        self.current_results_queue.append(result)

    def _add_fallback_result(self, smiles: str, identifier: str) -> None:
        """
        Add a fallback result when processing fails.

        Args:
            smiles: Original SMILES string to use as fallback
            identifier: Molecule identifier
        """
        assert isinstance(smiles, str)
        assert isinstance(identifier, str)
        assert len(smiles) > 0

        result = ProtonationResult(smiles=smiles, identifier=identifier, states="")
        self.current_results_queue.append(result)
        self.stats.fallback_used += 1

        logger.debug("Added fallback result for '{}'", smiles)

    def _format_statistics(self) -> str:
        """
        Format processing statistics for logging.

        Returns:
            Formatted statistics string
        """
        return (
            f"processed: {self.stats.molecules_processed}, "
            f"with_sites: {self.stats.molecules_with_sites}, "
            f"without_sites: {self.stats.molecules_without_sites}, "
            f"variants: {self.stats.total_variants_generated}, "
            f"validated: {self.stats.variants_validated}, "
            f"fallbacks: {self.stats.fallback_used}"
        )

    def stream_all(self) -> Generator[str, None, None]:
        """
        Stream all protonated SMILES as a generator.

        Yields:
            Protonated SMILES strings with metadata
        """
        try:
            while True:
                yield next(self)
        except StopIteration:
            return

    def to_list(self) -> list[str]:
        """
        Return all protonated SMILES as a list.

        Returns:
            List of protonated SMILES strings with metadata
        """
        return list(self.stream_all())

    def get_stats(self) -> dict[str, int | dict[str, int]]:
        """
        Get comprehensive processing statistics.

        Returns:
            Dictionary containing processing statistics
        """
        processor_stats = self.smiles_processor.get_stats()

        # Convert dataclass to dict for consistency
        protonation_stats = {
            "molecules_processed": self.stats.molecules_processed,
            "molecules_with_sites": self.stats.molecules_with_sites,
            "molecules_without_sites": self.stats.molecules_without_sites,
            "total_variants_generated": self.stats.total_variants_generated,
            "variants_validated": self.stats.variants_validated,
            "variants_rejected": self.stats.variants_rejected,
            "fallback_used": self.stats.fallback_used,
        }

        substructure_count = 0
        if self.pka_data and hasattr(self.pka_data, "_data"):
            substructure_count = len(self.pka_data._data)

        return {
            "processor": processor_stats,
            "protonation": protonation_stats,
            "total_substructures_loaded": substructure_count,
        }

    def reset_stats(self) -> None:
        """Reset all processing statistics to zero."""
        self.stats = ProtonationStats()
        logger.debug("Reset protonation statistics")


def protonate_smiles(
    smiles_input: str | Iterable[str] | Iterator[str],
    ph_min: float = 6.4,
    ph_max: float = 8.4,
    precision: float = 1.0,
    label_identifiers: bool = False,
    label_states: bool = False,
    max_variants: int = 128,
    validate_output: bool = True,
    **kwargs,
) -> list[str]:
    """
    Convenience function to protonate SMILES with explicit parameters.

    Args:
        smiles_input: SMILES string, file path, or iterable of SMILES
        ph_min: Minimum pH to consider
        ph_max: Maximum pH to consider
        precision: pKa precision factor
        include_states: When returning SMILES, format the string to include
            states.
        label_identifiers: When returning SMILES, format the string to
            include any identifier.
        max_variants: Maximum number of variants per input compound
        validate_output: Whether to validate generated SMILES
        **kwargs: Additional arguments for SMILESProcessor

    Returns:
        Generator of protonated SMILES strings
    """
    assert isinstance(ph_min, (int, float))
    assert isinstance(ph_max, (int, float))
    assert isinstance(precision, (int, float))
    assert isinstance(label_identifiers, bool)
    assert isinstance(label_states, bool)
    assert isinstance(max_variants, int)
    assert isinstance(validate_output, bool)

    protonator = Protonate(
        smiles_input=smiles_input,
        ph_min=float(ph_min),
        ph_max=float(ph_max),
        precision=float(precision),
        label_identifiers=label_identifiers,
        label_states=label_states,
        max_variants=max_variants,
        validate_output=validate_output,
        **kwargs,
    )

    return protonator.to_list()
