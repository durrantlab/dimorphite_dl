from typing import Generator

from dataclasses import dataclass, field
from enum import Enum

from dimorphite_dl.protonate.data import SubstructureDatum


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


@dataclass
class ProtonationSite:
    """Data class for protonation site information."""

    idx_atom: int
    """Atom index of protonation site"""
    substructure: SubstructureDatum = field(default_factory=SubstructureDatum)
    """Data of this protonation site matched by substructure"""

    def get_states(
        self, ph_min: float, ph_max: float, pka_stdev_prefactor: float
    ) -> Generator[ProtonationState]:
        for pka in self.substructure.pkas:
            yield pka.get_state(ph_min, ph_max, pka_stdev_prefactor)
