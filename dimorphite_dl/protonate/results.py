from dataclasses import dataclass

from dimorphite_dl.protonate.site import ProtonationState


@dataclass
class ProtonationResult:
    """Data class for protonation results."""

    smiles: str
    identifier: str
    states: ProtonationState | None = None

    def to_string(self, include_states: bool = False) -> str:
        """Convert to output string format."""
        if include_states and self.states:
            return f"{self.smiles}\t{self.identifier}\t{self.states.to_str()}"
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
