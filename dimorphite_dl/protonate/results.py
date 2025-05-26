from dataclasses import dataclass


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
