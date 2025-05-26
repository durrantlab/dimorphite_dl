from typing import Any

import importlib.resources as pkg_resources
from dataclasses import dataclass

from loguru import logger
from rdkit import Chem
from rdkit.Chem import Mol

from dimorphite_dl import data  # type: ignore
from dimorphite_dl.protonate.site import ProtonationState


@dataclass
class PkaData:
    site: int
    mean: float
    stdev: float
    state: ProtonationState | None = None

    def at_ph(self, ph_min: float, ph_max: float, pka_stdev_prefactor: float) -> None:
        """Computes the protonation state of at this pKa."""

        assert isinstance(ph_min, float)
        assert isinstance(ph_max, float)
        assert isinstance(pka_stdev_prefactor, float)

        stdev = pka_stdev_prefactor * self.stdev

        pka_min = self.mean - stdev
        pka_max = self.mean + stdev

        # This needs to be reassigned, and 'ERROR' should never make it past
        # the next set of checks.
        if (pka_min <= ph_max) and (ph_min <= pka_max):
            protonation_state = ProtonationState.BOTH
        elif self.mean > ph_max:
            protonation_state = ProtonationState.PROTONATED
        else:
            protonation_state = ProtonationState.DEPROTONATED

        self.state = protonation_state


@dataclass
class ProtonationSiteData:
    """Data class for protonation site data from database."""

    name: str
    """Name of the substructure"""
    smarts: str
    """SMARTS pattern for substructure"""
    pkas: list[PkaData] = []
    """pKa data observed from the dataset and possibly computed states"""
    mol: Mol | None = None
    """RDKit Mol"""

    def at_ph(self, ph_min: float, ph_max: float, pka_stdev_prefactor: float) -> None:
        for pka in self.pkas:
            pka.at_ph(ph_min, ph_max, pka_stdev_prefactor)


class SubstructureData:
    _data: list[ProtonationSiteData] = []
    """All loaded data for our protonation sites."""
    _data_at_ph: dict[tuple[float, float, float], list[ProtonationSiteData]] = {}
    """
    Cache of computed protonation sites based on `(ph_min, ph_max, ph_stdev)` keys.
    """

    @staticmethod
    def _load_lines() -> list[str]:
        """Load the substructure SMARTS file, filtering out comments and blank lines.

        Returns:
            List of valid SMARTS lines from the file.

        Raises:
            FileNotFoundError: If the substructure file cannot be found.
            IOError: If there are issues reading the file.
        """
        logger.trace("Loading substructure data from site_substructures.smarts")

        try:
            with pkg_resources.open_text(data, "site_substructures.smarts") as f:
                lines = []
                line_count = 0
                valid_count = 0

                for line in f:
                    line_count += 1
                    stripped = line.strip()

                    # Skip empty lines and comments
                    if stripped and not stripped.startswith("#"):
                        lines.append(stripped)
                        valid_count += 1

                logger.debug(
                    "Loaded {} valid SMARTS patterns from {} total lines",
                    valid_count,
                    line_count,
                )
                return lines

        except FileNotFoundError:
            logger.error("Could not find site_substructures.smarts file")
            raise
        except Exception as e:
            logger.error("Error reading substructure file: {}", str(e))
            raise IOError(f"Failed to read substructure file: {e}")

    @classmethod
    def _parse_substructure_line(cls, line: str) -> ProtonationSiteData:
        """Parse a single line from the substructure data file.

        Args:
            line: Line from the substructure file
            min_ph: Minimum pH
            max_ph: Maximum pH
            pka_stdev_prefactor: pKa standard range multiplier

        Returns:
            SubstructureData object.

        Notes:
            Below is an example line of the tab separated file.

            ```text
            *Azide	[N+0:1]=[N+:2]=[N+0:3]-[H]	2	4.65	0.07071067811865513
            ```

            This contains the following information separated by tabs.

            -   Name of the substructure. A `*` prefix indicates that it is an aromatic
                nitrogen that needs special treatment.
            -   [SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
                of this particular substructure.
            -   Data about the protonation site always in a set of threes. You can have
                more than one site.
                - The site index.
                - pKa mean
                - pKa standard deviation.
        """
        parts = line.split()
        if len(parts) < 3:
            logger.warning("Invalid line format (too few parts): '{}'", line)
            raise ValueError

        name = parts[0]
        logger.trace("Substructure name is {}", name)
        smarts = parts[1]
        logger.trace("Substructure SMARTS is {}", smarts)
        mol = cls._create_rdkit_mol(smarts)

        # Parse pKa ranges (groups of 3: site, mean, std)
        pka_data = cls._parse_pka_line(parts[2:])
        return ProtonationSiteData(name=name, smarts=smarts, pkas=pka_data, mol=mol)

    @classmethod
    def _create_rdkit_mol(cls, smarts: str) -> Chem.Mol:
        # Create mol object from SMARTS
        try:
            logger.trace("Attempting to make RDKit mol from SMARTS")
            mol = Chem.MolFromSmarts(smarts)
            if mol is None:
                logger.warning("Invalid SMARTS pattern: {}", smarts)
                raise ValueError
        except Exception as e:
            logger.warning("Error creating mol from SMARTS '{}' : {}", smarts, str(e))
            raise ValueError
        return mol

    @classmethod
    def _parse_pka_line(cls, line_parts: list[str]) -> list[PkaData]:
        if len(line_parts) % 3 != 0:
            logger.warning(
                "Invalid pKa data format, expected groups of 3, got {}", len(line_parts)
            )
            raise ValueError

        pka_data = []
        for i in range(0, len(line_parts), 3):
            try:
                site = int(line_parts[i])
                mean = float(line_parts[i + 1])
                stdev = float(line_parts[i + 2])
                pka_data.append(PkaData(site=site, mean=mean, stdev=stdev))
            except (ValueError, IndexError) as e:
                logger.warning("Error parsing pKa data: {}", line_parts)
                raise ValueError from e
        return pka_data

    def get_protonation_states(
        self, ph_min: float, ph_max: float, pka_stdev_prefactor: float
    ) -> list[ProtonationSiteData]:
        """Parse a single line from the substructure file.

        Args:
            line: Line from the substructure file
            ph_max: Minimum pH
            ph_max: Maximum pH
            pka_stdev_prefactor: pKa standard range multiplier

        Returns:
            SubstructureData object or None if parsing fails
        """
        key_cache = (ph_min, ph_max, pka_stdev_prefactor)
        if key_cache in self._data_at_ph.keys():
            return self._data_at_ph[key_cache]
        else:
            pka_at_ph = self.compute_prot_states(ph_min, ph_max, pka_stdev_prefactor)
            self._data_at_ph[key_cache] = pka_at_ph
            return pka_at_ph

    def compute_prot_states(
        self, ph_min: float, ph_max: float, pka_stdev_prefactor: float
    ) -> list[ProtonationSiteData]:
        prot_states = []
        for prot_site in self._data:
            std = std_base * pka_stdev_prefactor
            protonation_state = define_protonation_state(mean, std, min_ph, max_ph)
            logger.trace("Calculated protonation state of {}", protonation_state)
            prot_states.append([site, protonation_state])

        return SubstructureData(
            name=name, smart=smart, mol=mol, prot_states_for_pH=prot_states
        )


def load_protonation_substructs_calc_state_for_ph(
    min_ph: float = 6.4, max_ph: float = 8.4, pka_stdev_prefactor: float = 1.0
) -> list[SubstructureData]:
    """Load protonation substructures with calculated states for pH range.

    Args:
        min_ph: Lower bound of pH range (default: 6.4)
        max_ph: Upper bound of pH range (default: 8.4)
        pka_stdev_prefactor: Precision factor for pKa calculations (default: 1.0)

    Returns:
        List of SubstructureData objects with protonation information.

    Raises:
        ValueError: If pH range is invalid.
    """
    if min_ph > max_ph:
        raise ValueError(f"Invalid pH range: min_ph ({min_ph}) >= max_ph ({max_ph})")

    if pka_stdev_prefactor < 0:
        raise ValueError(
            f"pka_stdev_prefactor must be positive, got: {pka_stdev_prefactor}"
        )

    logger.info(
        "Loading protonation substructures for pH range {:.1f}-{:.1f}", min_ph, max_ph
    )

    subs = []
    failed_count = 0

    for line_num, line in enumerate(load_substructure_smarts_file(), 1):
        try:
            sub_data = _parse_substructure_line(
                line, min_ph, max_ph, pka_stdev_prefactor
            )
            if sub_data:
                subs.append(sub_data)
            else:
                failed_count += 1
        except Exception as e:
            logger.warning("Failed to parse line {}: '{}' - {}", line_num, line, str(e))
            failed_count += 1

    logger.info(
        "Successfully loaded {} substructures, {} failed", len(subs), failed_count
    )
    return subs


def validate_substructure_data(subs: list[SubstructureData]) -> dict[str, Any]:
    """Validate loaded substructure data and return statistics.

    Args:
        subs: List of substructure data

    Returns:
        Dictionary with validation statistics
    """
    stats: dict[str, int | dict[str, int]] = {
        "total_substructures": len(subs),
        "valid_mols": 0,
        "invalid_mols": 0,
        "total_sites": 0,
        "state_distribution": {"PROTONATED": 0, "DEPROTONATED": 0, "BOTH": 0},
    }

    for sub in subs:
        if sub.mol is not None:
            stats["valid_mols"] += 1
        else:
            stats["invalid_mols"] += 1

        stats["total_sites"] += len(sub.prot_states_for_pH)

        for site_info in sub.prot_states_for_pH:
            state = site_info[1]
            if state in stats["state_distribution"]:
                stats["state_distribution"][state] += 1

    return stats
