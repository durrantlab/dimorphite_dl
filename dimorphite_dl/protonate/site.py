from dataclasses import dataclass
from enum import Enum


class ProtonationState(Enum):
    DEPROTONATED = 0
    PROTONATED = 1
    BOTH = 2


@dataclass
class ProtonationSite:
    """Data class for protonation site information."""

    atom_idx: int
    target_state: ProtonationState
    site_name: str
