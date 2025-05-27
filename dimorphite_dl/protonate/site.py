from typing import Generator

from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem


class ProtonationState(Enum):
    UNKNOWN = 0
    DEPROTONATED = 1
    PROTONATED = 2
    BOTH = 3

    def to_str(self) -> str:
        if self == ProtonationState.DEPROTONATED:
            return "DEPROTONATED"
        elif self == ProtonationState.PROTONATED:
            return "PROTONATED"
        elif self == ProtonationState.BOTH:
            return "BOTH"
        else:
            return "UNKNOWN"

    def get_charges(self) -> list[int]:
        if self == ProtonationState.DEPROTONATED:
            return [-1]
        elif self == ProtonationState.PROTONATED:
            return [0]
        elif self == ProtonationState.BOTH:
            return [-1, 0]
        else:
            return []


class PKaDatum:
    def __init__(self, site, mean, stdev):
        self.site = site
        self.mean = mean
        self.stdev = stdev

    def get_state(
        self, ph_min: float, ph_max: float, precision: float
    ) -> ProtonationState:
        """Computes the protonation state at this pH."""

        assert isinstance(ph_min, float)
        assert isinstance(ph_max, float)
        assert isinstance(precision, float)

        stdev = precision * self.stdev

        pka_min = self.mean - stdev
        pka_max = self.mean + stdev

        # This needs to be reassigned, and 'ERROR' should never make it past
        # the next set of checks.
        if (pka_min <= ph_max) and (ph_min <= pka_max):
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
    """Data class for protonation site data from database."""

    name: str = ""
    """Name of the substructure"""
    smarts: str = ""
    """SMARTS pattern for substructure"""
    mol: Chem.Mol | None = None
    """RDKit Mol of substructure"""
    pkas: list[PKaDatum] = field(default_factory=list)
    """pKa data observed from the dataset and possibly computed states"""

    def at_ph(self, ph_min: float, ph_max: float, precision: float) -> None:
        for pka in self.pkas:
            pka.at_ph(ph_min, ph_max, precision)


@dataclass
class ProtonationSite:
    """Data class for protonation site information."""

    idx_atom: int
    """Atom index of protonation site"""
    substructure: SubstructureDatum = field(default_factory=SubstructureDatum)
    """Data of this protonation site matched by substructure"""

    def get_states(
        self, ph_min: float, ph_max: float, precision: float
    ) -> Generator[ProtonationState]:
        for pka in self.substructure.pkas:
            yield pka.get_state(ph_min, ph_max, precision)
