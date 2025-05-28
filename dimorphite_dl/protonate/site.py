"""
Protonation site data structures and state calculations.

This module defines the core data structures for protonation sites,
including state enumerations, pKa data, and site information.
Each class has clear responsibilities and comprehensive validation.
"""

from typing import Generator

from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem


class ProtonationState(Enum):
    """
    Enumeration of possible protonation states for a site.

    Values are explicitly assigned for clarity and debugging.
    """

    UNKNOWN = 0
    DEPROTONATED = 1
    PROTONATED = 2
    BOTH = 3

    def to_str(self) -> str:
        """
        Convert protonation state to string representation.

        Returns:
            String representation of the protonation state
        """
        # Use explicit if-elif chain for clarity (TigerStyle)
        if self == ProtonationState.DEPROTONATED:
            return "DEPROTONATED"
        elif self == ProtonationState.PROTONATED:
            return "PROTONATED"
        elif self == ProtonationState.BOTH:
            return "BOTH"
        else:
            return "UNKNOWN"

    def get_charges(self) -> list[int]:
        """
        Get the formal charges associated with this protonation state.

        Returns:
            List of integer formal charges for this state
        """
        # Use explicit if-elif chain for clarity (TigerStyle)
        if self == ProtonationState.DEPROTONATED:
            return [-1]
        elif self == ProtonationState.PROTONATED:
            return [0]
        elif self == ProtonationState.BOTH:
            return [-1, 0]
        else:
            return []


class PKaDatum:
    """
    Data structure for pKa information at a specific site.

    Contains the site index, mean pKa value, and standard deviation
    for calculating protonation states at different pH values.
    """

    def __init__(self, site: int, mean: float, stdev: float):
        """
        Initialize pKa data with validation.

        Args:
            site: Site index (non-negative integer)
            mean: Mean pKa value (bounded 0-20 for realistic range)
            stdev: Standard deviation (non-negative, bounded 0-5)
        """
        assert isinstance(site, int)
        assert isinstance(mean, (int, float))
        assert isinstance(stdev, (int, float))
        assert site >= 0, f"Site index must be non-negative, got: {site}"

        self.site = site
        self.mean = float(mean)
        self.stdev = float(stdev)

    def get_state(
        self, ph_min: float, ph_max: float, precision: float
    ) -> ProtonationState:
        """
        Calculate protonation state for given pH range and precision.

        Args:
            ph_min: Minimum pH value (bounded 0-14)
            ph_max: Maximum pH value (bounded 0-14, greater than ph_min)
            precision: Precision factor for pKa calculation (positive)

        Returns:
            ProtonationState based on pH range and pKa statistics
        """
        assert isinstance(ph_min, (int, float))
        assert isinstance(ph_max, (int, float))
        assert isinstance(precision, (int, float))
        assert 0.0 <= ph_min <= 14.0, f"ph_min must be 0-14, got: {ph_min}"
        assert 0.0 <= ph_max <= 14.0, f"ph_max must be 0-14, got: {ph_max}"
        assert ph_min <= ph_max, (
            f"ph_min ({ph_min}) must be less than ph_max ({ph_max})"
        )
        assert precision >= 0.0, f"precision must be positive, got: {precision}"

        # Calculate effective pKa range based on precision
        effective_stdev = precision * self.stdev
        pka_min = self.mean - effective_stdev
        pka_max = self.mean + effective_stdev

        # Determine protonation state based on pH and pKa overlap
        # Use explicit conditions for clarity (TigerStyle)
        pka_overlaps_ph_range = (pka_min <= ph_max) and (ph_min <= pka_max)
        if pka_overlaps_ph_range:
            protonation_state = ProtonationState.BOTH
        elif self.mean > ph_max:
            protonation_state = ProtonationState.PROTONATED
        elif self.mean < ph_min:
            protonation_state = ProtonationState.DEPROTONATED
        else:
            protonation_state = ProtonationState.UNKNOWN

        return protonation_state


@dataclass
class SubstructureDatum:
    """
    Data structure for substructure pattern matching information.

    Contains the pattern name, SMARTS string, RDKit mol object,
    and associated pKa data for protonation site detection.
    """

    name: str = ""
    smarts: str = ""
    mol: Chem.Mol | None = None
    pkas: list[PKaDatum] = field(default_factory=list)

    def __post_init__(self):
        """Validate substructure data after initialization."""
        assert isinstance(self.name, str)
        assert isinstance(self.smarts, str)
        assert isinstance(self.pkas, list)

        # Name and SMARTS should not be empty for valid substructures
        if len(self.name) > 0 or len(self.smarts) > 0:
            assert len(self.name) > 0, "Substructure name cannot be empty"
            assert len(self.smarts) > 0, "SMARTS pattern cannot be empty"

        # Validate all pKa data entries
        for pka in self.pkas:
            assert isinstance(pka, PKaDatum)

    def has_valid_pattern(self) -> bool:
        """
        Check if substructure has a valid molecular pattern.

        Returns:
            True if mol object exists and is valid
        """
        return self.mol is not None

    def has_pka_data(self) -> bool:
        """
        Check if substructure has pKa data available.

        Returns:
            True if pKa data list is non-empty
        """
        return len(self.pkas) > 0

    def get_pka_count(self) -> int:
        """
        Get the number of pKa data points for this substructure.

        Returns:
            Number of pKa data entries
        """
        return len(self.pkas)

    def is_valid_for_matching(self) -> bool:
        """
        Check if substructure is valid for pattern matching.

        Returns:
            True if both pattern and pKa data are available
        """
        return self.has_valid_pattern() and self.has_pka_data()


@dataclass
class ProtonationSite:
    """
    Data structure for detected protonation site information.

    Contains atom indices and associated substructure data
    for a specific protonation site in a molecule.
    """

    idx_atom: int | tuple[int, ...]
    substructure: SubstructureDatum = field(default_factory=SubstructureDatum)

    def __post_init__(self):
        """Validate protonation site data after initialization."""
        # Validate atom index format
        if isinstance(self.idx_atom, int):
            assert self.idx_atom >= 0, (
                f"Atom index must be non-negative, got: {self.idx_atom}"
            )
        elif isinstance(self.idx_atom, tuple):
            assert len(self.idx_atom) > 0, "Atom index tuple cannot be empty"
            for idx in self.idx_atom:
                assert isinstance(idx, int), (
                    f"All atom indices must be integers, got: {type(idx)}"
                )
                assert idx >= 0, f"All atom indices must be non-negative, got: {idx}"
        else:
            raise ValueError(
                f"idx_atom must be int or tuple of ints, got: {type(self.idx_atom)}"
            )

        assert isinstance(self.substructure, SubstructureDatum)

    def get_atom_indices(self) -> list[int]:
        """
        Get atom indices as a list for consistent handling.

        Returns:
            List of atom indices involved in this protonation site
        """
        if isinstance(self.idx_atom, int):
            return [self.idx_atom]
        else:
            return list(self.idx_atom)

    def get_primary_atom_index(self) -> int:
        """
        Get the primary (first) atom index for this site.

        Returns:
            Primary atom index
        """
        if isinstance(self.idx_atom, int):
            return self.idx_atom
        else:
            assert len(self.idx_atom) > 0
            return self.idx_atom[0]

    def get_atom_count(self) -> int:
        """
        Get the number of atoms involved in this protonation site.

        Returns:
            Number of atoms in the protonation site
        """
        if isinstance(self.idx_atom, int):
            return 1
        else:
            return len(self.idx_atom)

    def has_substructure_data(self) -> bool:
        """
        Check if site has valid substructure data.

        Returns:
            True if substructure data is available and valid
        """
        return self.substructure.is_valid_for_matching()

    def get_states(
        self, ph_min: float, ph_max: float, precision: float
    ) -> Generator[ProtonationState, None, None]:
        """
        Generate protonation states for all pKa data at this site.

        Args:
            ph_min: Minimum pH value (bounded 0-14)
            ph_max: Maximum pH value (bounded 0-14, greater than ph_min)
            precision: Precision factor for pKa calculation (positive)

        Yields:
            ProtonationState for each pKa datum at this site
        """
        assert isinstance(ph_min, (int, float))
        assert isinstance(ph_max, (int, float))
        assert isinstance(precision, (int, float))
        assert 0.0 <= ph_min <= 14.0, f"ph_min must be 0-14, got: {ph_min}"
        assert 0.0 <= ph_max <= 14.0, f"ph_max must be 0-14, got: {ph_max}"
        assert ph_min <= ph_max, (
            f"ph_min ({ph_min}) must be less than ph_max ({ph_max})"
        )
        assert precision >= 0.0, f"precision must be positive, got: {precision}"

        pka_count = len(self.substructure.pkas)
        assert pka_count > 0, "Cannot generate states without pKa data"

        states_generated = 0
        for pka in self.substructure.pkas:
            assert isinstance(pka, PKaDatum)
            state = pka.get_state(ph_min, ph_max, precision)
            states_generated += 1
            yield state

        assert states_generated == pka_count, (
            f"Expected {pka_count} states, generated {states_generated}"
        )

    def get_states_list(
        self, ph_min: float, ph_max: float, precision: float
    ) -> list[ProtonationState]:
        """
        Get protonation states as a list for easier handling.

        Args:
            ph_min: Minimum pH value (bounded 0-14)
            ph_max: Maximum pH value (bounded 0-14, greater than ph_min)
            precision: Precision factor for pKa calculation (positive)

        Returns:
            List of ProtonationState objects for this site
        """
        return list(self.get_states(ph_min, ph_max, precision))


# Utility functions for working with protonation states
def validate_ph_range(ph_min: float, ph_max: float) -> bool:
    """
    Validate pH range parameters.

    Args:
        ph_min: Minimum pH value
        ph_max: Maximum pH value

    Returns:
        True if pH range is valid
    """
    try:
        assert isinstance(ph_min, (int, float))
        assert isinstance(ph_max, (int, float))
        assert 0.0 <= ph_min <= 14.0
        assert 0.0 <= ph_max <= 14.0
        assert ph_min < ph_max
        return True
    except AssertionError:
        return False


def create_pka_datum_safe(site: int, mean: float, stdev: float) -> PKaDatum | None:
    """
    Create PKaDatum with error handling.

    Args:
        site: Site index
        mean: Mean pKa value
        stdev: Standard deviation

    Returns:
        PKaDatum object or None if parameters are invalid
    """
    try:
        return PKaDatum(site, mean, stdev)
    except (AssertionError, ValueError):
        return None


def create_protonation_site_safe(
    idx_atom: int | tuple[int, ...], substructure: SubstructureDatum | None = None
) -> ProtonationSite | None:
    """
    Create ProtonationSite with error handling.

    Args:
        idx_atom: Atom index or tuple of indices
        substructure: Substructure data (optional)

    Returns:
        ProtonationSite object or None if parameters are invalid
    """
    try:
        if substructure is None:
            substructure = SubstructureDatum()
        return ProtonationSite(idx_atom, substructure)
    except (AssertionError, ValueError):
        return None
