from typing import Generator, Iterable, Iterator

from dataclasses import dataclass
from pathlib import Path

from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem

from dimorphite_dl.io import SMILESProcessor, SMILESRecord
from dimorphite_dl.mol import smiles_to_mol
from dimorphite_dl.protonate import PKaData
from dimorphite_dl.protonate.detect import ProtonationSiteDetector
from dimorphite_dl.protonate.substruct import protonate_site


@dataclass
class ProtonationResult:
    """Data class for protonation results."""

    smiles: str
    identifier: str
    states: str = ""

    def to_string(self, include_states: bool = False) -> str:
        """Convert to output string format."""
        if include_states and self.states:
            return f"{self.smiles}\t{self.identifier}\t{self.states}"
        return f"{self.smiles}\t{self.identifier}"


@dataclass
class ProtonationStats:
    """Statistics for protonation processing."""

    molecules_processed: int = 0
    total_variants_generated: int = 0
    variants_validated: int = 0
    variants_rejected: int = 0
    molecules_with_sites: int = 0
    molecules_without_sites: int = 0
    fallback_used: int = 0


class Protonate:
    """A generator class for protonating SMILES strings, one at a time."""

    def __init__(
        self,
        smiles_input: str | Iterable[str] | Iterator[str] | None = None,
        ph_min: float = 6.4,
        ph_max: float = 8.4,
        pka_stdev_prefactor: float = 1.0,
        label_states: bool = False,
        max_variants: int = 128,
        validate_output: bool = True,
        **smiles_processor_kwargs,
    ):
        """Initialize the generator.

        Args:
            smiles_input: SMILES string, file path, or iterable of SMILES
            ph_min: Minimum pH to consider (default: 6.4)
            ph_max: Maximum pH to consider (default: 8.4)
            pka_stdev_prefactor: pKa precision factor (default: 1.0)
            label_states: Whether to label protonated SMILES with target state
            max_variants: Maximum number of variants per input compound
            validate_output: Whether to validate generated SMILES (default: True)
            **smiles_processor_kwargs: Additional arguments for SMILESProcessor
        """
        # Validate input parameters
        if ph_min > ph_max:
            raise ValueError(f"ph_min ({ph_min}) must be less than ph_max ({ph_max})")
        if pka_stdev_prefactor < 0:
            raise ValueError(
                f"pka_stdev_prefactor must be positive, got: {pka_stdev_prefactor}"
            )
        if max_variants <= 0:
            raise ValueError(f"max_variants must be positive, got: {max_variants}")

        self.smiles_input = smiles_input
        self.ph_min = ph_min
        self.ph_max = ph_max
        self.pka_stdev_prefactor = pka_stdev_prefactor
        self.label_states = label_states
        self.max_variants = max_variants
        self.validate_output = validate_output

        logger.info(
            "Initializing with pH range {:.1f} to {:.1f}, precision: {:.1f}",
            ph_min,
            ph_max,
            pka_stdev_prefactor,
        )

        self.smiles_processor = SMILESProcessor(**smiles_processor_kwargs)
        self.stats = ProtonationStats()

        # Queue to store the protonated SMILES strings for current molecule
        self.cur_prot_results: list[ProtonationResult] = []

        self.pka_data = PKaData()

        # Initialize the SMILES stream
        self._smiles_stream: Iterator[SMILESRecord] | None = None
        self._initialize_stream()

    def _initialize_stream(self) -> None:
        """Initialize the SMILES stream from the input."""
        if self.smiles_input is not None:
            try:
                self._smiles_stream = iter(
                    self.smiles_processor.stream(self.smiles_input)
                )
                logger.debug("Initialized SMILES stream")
            except Exception as e:
                logger.error("Failed to initialize SMILES stream: {}", str(e))
                raise ValueError(f"Failed to initialize SMILES stream: {e}")
        else:
            raise ValueError("No SMILES input provided")

    def __iter__(self):
        """Returns this generator object."""
        return self

    def __next__(self) -> str:
        """Return the next protonated SMILES string.

        Returns:
            A string containing the protonated SMILES with metadata.

        Raises:
            StopIteration: When no more SMILES are available to process.
        """
        # If there are any results in the queue, return the first one
        if self.cur_prot_results:
            result = self.cur_prot_results.pop(0)
            return result.to_string(include_states=self.label_states)

        # If no current results, process the next input SMILES
        if self._smiles_stream is None:
            raise StopIteration("No SMILES stream available")

        try:
            smiles_record = next(self._smiles_stream)
            self._process_smiles_record(smiles_record)

            # After processing, there should be results in the queue
            if self.cur_prot_results:
                result = self.cur_prot_results.pop(0)
                return result.to_string(include_states=self.label_states)
            else:
                # If processing didn't generate any results, try next record
                return self.__next__()

        except StopIteration:
            # No more input SMILES
            logger.info("Processing complete. Final stats: {}", self._format_stats())
            raise StopIteration("No more SMILES to process")

    def _process_smiles_record(self, smiles_record: SMILESRecord) -> None:
        """Process a single SMILES record and generate protonated variants.

        Args:
            smiles_record: SMILESRecord object containing SMILES and metadata
        """
        orig_smi = smiles_record.smiles
        identifier = smiles_record.identifier or ""

        orig_smi = self._neutralize_smiles(orig_smi)

        self.stats.molecules_processed += 1

        if orig_smi is None:
            return

        # Track valid SMILES for fallback purposes
        valid_smiles_history = [orig_smi]

        try:
            mol, sites = ProtonationSiteDetector.find(orig_smi)

            if (mol is None) or (len(sites) == 0):
                logger.warning("Could not find protonation sites for {}", orig_smi)
                self._add_fallback_result(orig_smi, identifier)
                return

            self.stats.molecules_with_sites += 1
            new_mols = [mol]

            # Process each protonation site
            for i, site in enumerate(sites):
                try:
                    new_mols = protonate_site(
                        new_mols,
                        site,
                        self.ph_min,
                        self.ph_max,
                        self.pka_stdev_prefactor,
                    )
                except Exception as e:
                    logger.warning(
                        "Error processing site {} for '{}': {}", i + 1, orig_smi, str(e)
                    )
                    continue

                # Limit variants if necessary
                if len(new_mols) > self.max_variants:
                    new_mols = new_mols[: self.max_variants]
                    logger.warning(
                        "Limited variants to {} for '{}' at site {}",
                        self.max_variants,
                        orig_smi,
                        i + 1,
                    )

                # Add valid SMILES to history for potential fallback
                try:
                    valid_smiles_history.extend(
                        [Chem.MolToSmiles(m) for m in new_mols if m is not None]
                    )
                except Exception as e:
                    logger.debug("Error generating SMILES for history: {}", str(e))

            else:
                logger.debug("No protonation sites found for '{}'", orig_smi)
                self.stats.molecules_without_sites += 1

                # Remove hydrogens since protonate_site was never called
                try:
                    mol_used_to_idx_sites = Chem.RemoveHs(mol)
                    if mol_used_to_idx_sites is not None:
                        new_mols = [mol_used_to_idx_sites]
                        valid_smiles_history.append(
                            Chem.MolToSmiles(mol_used_to_idx_sites)
                        )
                except Exception as e:
                    logger.warning(
                        "Error removing hydrogens from '{}': {}", orig_smi, str(e)
                    )

            # Generate canonical SMILES and remove duplicates
            canonical_smiles = self._generate_canonical_smiles(new_mols, orig_smi)

            # Validate generated SMILES if requested
            if self.validate_output:
                validated_smiles = self._validate_smiles_list(
                    canonical_smiles, orig_smi
                )
            else:
                validated_smiles = canonical_smiles

            # If no valid SMILES remain, use fallback
            if not validated_smiles:
                logger.warning(
                    "No valid SMILES generated for '{}', using fallback", orig_smi
                )
                validated_smiles = self._get_fallback_smiles(valid_smiles_history)
                self.stats.fallback_used += 1

            # Create results with state information
            self._create_protonation_results(validated_smiles, identifier, sites)

        except Exception as e:
            logger.error("Error processing SMILES '{}': {}", orig_smi, str(e))
            self._add_fallback_result(orig_smi, identifier)

    def _neutralize_smiles(self, smi: str) -> str | None:
        """All molecules should be neuralized to the extent possible. The user
        should not be allowed to specify the valence of the atoms in most cases.

        Args:
            mol: The rdkit Mol objet to be neutralized.

        Returns:
            The neutralized Mol object.
        """
        logger.debug("Neutralizing {}", smi)
        # Get the reaction data
        rxn_data = [
            [
                "[Ov1-1:1]",
                "[Ov2+0:1]-[H]",
            ],  # To handle O- bonded to only one atom (add hydrogen).
            [
                "[#7v4+1:1]-[H]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to a hydrogen (remove hydrogen).
            [
                "[Ov2-:1]",
                "[Ov2+0:1]",
            ],  # To handle O- bonded to two atoms. Should not be Negative.
            [
                "[#7v3+1:1]",
                "[#7v3+0:1]",
            ],  # To handle N+ bonded to three atoms. Should not be positive.
            [
                "[#7v2-1:1]",
                "[#7+0:1]-[H]",
            ],  # To handle N- Bonded to two atoms. Add hydrogen.
            # ['[N:1]=[N+0:2]=[N:3]-[H]', '[N:1]=[N+1:2]=[N+0:3]-[H]'],  # To handle bad azide. Must be
            # protonated. (Now handled
            # elsewhere, before SMILES
            # converted to Mol object.)
            [
                "[H]-[N:1]-[N:2]#[N:3]",
                "[N:1]=[N+1:2]=[N:3]-[H]",
            ],  # To handle bad azide. R-N-N#N should
            # be R-N=[N+]=N
        ]

        # Add substructures and reactions (initially none)
        for i, rxn_datum in enumerate(rxn_data):
            rxn_data[i].append(Chem.MolFromSmarts(rxn_datum[0]))
            rxn_data[i].append(None)

        mol = smiles_to_mol(smi)
        # Add hydrogens (respects valence, so incomplete).
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol)

        while True:  # Keep going until all these issues have been resolved.
            current_rxn = None  # The reaction to perform.
            current_rxn_str = None

            for i, rxn_datum in enumerate(rxn_data):
                (
                    reactant_smarts,
                    product_smarts,
                    substruct_match_mol,
                    rxn_placeholder,
                ) = rxn_datum
                if mol.HasSubstructMatch(substruct_match_mol):
                    if rxn_placeholder is None:
                        current_rxn_str = reactant_smarts + ">>" + product_smarts
                        logger.trace("Defining reaction: {}", current_rxn_str)
                        current_rxn = AllChem.ReactionFromSmarts(current_rxn_str)
                        rxn_data[i][3] = current_rxn  # Update the placeholder.
                    else:
                        current_rxn = rxn_data[i][3]
                    break

            # Perform the reaction if necessary
            if current_rxn is None:  # No reaction left, so break out of while loop.
                break
            else:
                mol = current_rxn.RunReactants((mol,))[0][0]
                mol.UpdatePropertyCache(strict=False)  # Update valences

        # The mols have been altered from the reactions described above, we
        # need to resanitize them. Make sure aromatic rings are shown as such
        # This catches all RDKit Errors. without the catchError and
        # sanitizeOps the Chem.SanitizeMol can crash the program.
        sanitize_string = Chem.SanitizeMol(
            mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors=True
        )
        if sanitize_string.name == "SANITIZE_NONE":
            return Chem.MolToSmiles(mol)
        else:
            return None

    def _generate_canonical_smiles(
        self, mols: list[Chem.Mol], orig_smi: str
    ) -> list[str]:
        """Generate canonical SMILES from molecule objects."""
        try:
            smiles_set = set()
            for mol in mols:
                if mol is not None:
                    try:
                        canonical = Chem.MolToSmiles(
                            mol, isomericSmiles=True, canonical=True
                        )
                        if canonical:
                            smiles_set.add(canonical)
                    except Exception as e:
                        logger.debug("Error generating canonical SMILES: {}", str(e))
                        continue

            smiles_list = list(smiles_set)
            self.stats.total_variants_generated += len(smiles_list)
            logger.debug(
                "Generated {} unique canonical SMILES for '{}'",
                len(smiles_list),
                orig_smi,
            )
            return smiles_list

        except Exception as e:
            logger.warning(
                "Error in canonical SMILES generation for '{}': {}", orig_smi, str(e)
            )
            return []

    def _validate_smiles_list(self, smiles_list: list[str], orig_smi: str) -> list[str]:
        """Validate a list of SMILES strings."""
        if not smiles_list:
            return []

        logger.debug("Validating {} SMILES for '{}'", len(smiles_list), orig_smi)
        valid_smiles = []

        for smi in smiles_list:
            if smiles_to_mol(smi) is not None:
                valid_smiles.append(smi)
                self.stats.variants_validated += 1
            else:
                self.stats.variants_rejected += 1
                logger.debug("Rejected invalid SMILES: '{}'", smi)

        logger.debug(
            "Validated {}/{} SMILES for '{}'",
            len(valid_smiles),
            len(smiles_list),
            orig_smi,
        )
        return valid_smiles

    def _get_fallback_smiles(self, valid_history: list[str]) -> list[str]:
        """Get fallback SMILES from the history of valid SMILES."""
        # Try to find a valid SMILES from the history, starting from most recent
        for smi in reversed(valid_history):
            if smiles_to_mol(smi) is not None:
                logger.debug("Using fallback SMILES: '{}'", smi)
                return [smi]

        # If all else fails, return the original
        logger.warning("All fallback options failed, using original SMILES")
        return [valid_history[0]] if valid_history else []

    def _create_protonation_results(
        self, smiles_list: list[str], identifier: str, sites: list[int]
    ) -> None:
        """Create ProtonationResult objects from SMILES list."""
        if not smiles_list:
            return

        # Generate states string if needed
        states_str = ""
        if self.label_states and sites:
            states_str = "\t".join([site.target_state for site in sites])

        # Create results
        for smi in smiles_list:
            result = ProtonationResult(
                smiles=smi, identifier=identifier, states=states_str
            )
            self.cur_prot_results.append(result)

    def _add_fallback_result(self, smiles: str, identifier: str) -> None:
        """Add a fallback result when processing fails."""
        result = ProtonationResult(smiles=smiles, identifier=identifier)
        self.cur_prot_results.append(result)
        self.stats.fallback_used += 1

    def _format_stats(self) -> str:
        """Format statistics for logging."""
        return (
            f"processed: {self.stats.molecules_processed}, "
            f"with_sites: {self.stats.molecules_with_sites}, "
            f"variants: {self.stats.total_variants_generated}, "
            f"validated: {self.stats.variants_validated}, "
            f"fallbacks: {self.stats.fallback_used}"
        )

    def stream_all(self) -> Generator[str, None, None]:
        """Stream all protonated SMILES as a generator.

        Yields:
            Protonated SMILES strings with metadata
        """
        try:
            while True:
                yield next(self)
        except StopIteration:
            return

    def to_list(self) -> list[str]:
        """Return all protonated SMILES as a list.

        Returns:
            List of protonated SMILES strings with metadata
        """
        return list(self.stream_all())

    def get_stats(self) -> dict[str, int | dict[str, int]]:
        """Get comprehensive processing statistics.

        Returns:
            Dictionary containing processing statistics
        """
        processor_stats = self.smiles_processor.get_stats()

        protonation_stats = {
            "molecules_processed": self.stats.molecules_processed,
            "molecules_with_sites": self.stats.molecules_with_sites,
            "molecules_without_sites": self.stats.molecules_without_sites,
            "total_variants_generated": self.stats.total_variants_generated,
            "variants_validated": self.stats.variants_validated,
            "variants_rejected": self.stats.variants_rejected,
            "fallback_used": self.stats.fallback_used,
        }

        return {
            "processor": processor_stats,
            "protonation": protonation_stats,
            "total_substructures_loaded": len(self.pka_data._data)
            if self.pka_data
            else 0,
        }

    def reset_stats(self) -> None:
        """Reset processing statistics."""
        self.stats = ProtonationStats()
        logger.debug("Reset protonation statistics")


# Convenience functions for backward compatibility and ease of use
def protonate_smiles(
    smiles_input: str | Iterable[str] | Iterator[str] | Path,
    ph_min: float = 6.4,
    ph_max: float = 8.4,
    pka_stdev_prefactor: float = 1.0,
    label_states: bool = False,
    max_variants: int = 128,
    validate_output: bool = True,
    **kwargs,
) -> Generator[str, None, None]:
    """Convenience function to protonate SMILES.

    Args:
        smiles_input: SMILES string, file path, or iterable of SMILES
        ph_min: Minimum pH to consider
        ph_max: Maximum pH to consider
        pka_stdev_prefactor: pKa precision factor
        label_states: Whether to label protonated SMILES with target state
        max_variants: Maximum number of variants per input compound
        validate_output: Whether to validate generated SMILES
        **kwargs: Additional arguments for SMILESProcessor

    Returns:
        Generator of protonated SMILES strings
    """
    protonator = Protonate(
        smiles_input=smiles_input,
        ph_min=ph_min,
        ph_max=ph_max,
        pka_stdev_prefactor=pka_stdev_prefactor,
        label_states=label_states,
        max_variants=max_variants,
        validate_output=validate_output,
        **kwargs,
    )

    return protonator.stream_all()


def protonate_smiles_batch(
    smiles_list: list[str],
    ph_min: float = 6.4,
    ph_max: float = 8.4,
    pka_stdev_prefactor: float = 1.0,
    label_states: bool = False,
    max_variants: int = 128,
    validate_output: bool = True,
) -> dict[str, list[str] | dict[str, int]]:
    """Protonate a batch of SMILES and return results with statistics.

    Args:
        smiles_list: List of SMILES strings to protonate
        ph_min: Minimum pH to consider
        ph_max: Maximum pH to consider
        pka_stdev_prefactor: pKa precision factor
        label_states: Whether to label protonated SMILES with target state
        max_variants: Maximum number of variants per input compound
        validate_output: Whether to validate generated SMILES

    Returns:
        Dictionary with 'results' (list of protonated SMILES) and 'stats'
    """
    results = []
    protonator = Protonate(
        smiles_input=smiles_list,
        ph_min=ph_min,
        ph_max=ph_max,
        pka_stdev_prefactor=pka_stdev_prefactor,
        label_states=label_states,
        max_variants=max_variants,
        validate_output=validate_output,
    )

    for result in protonator.stream_all():
        results.append(result)

    return {"results": results, "stats": protonator.get_stats()}
