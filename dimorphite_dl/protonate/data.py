import importlib.resources as pkg_resources
from collections.abc import Iterator

from loguru import logger
from rdkit import Chem

from dimorphite_dl.protonate.site import PKaDatum, SubstructureDatum


class PKaData:
    _data: list[SubstructureDatum] = []
    """All loaded data for our protonation substructures."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._load_data()
        return cls._instance

    @classmethod
    def _load_data(cls) -> None:
        lines = cls._load_lines()
        data = []
        for line in lines:
            data.append(cls._parse_substructure_line(line))
        cls._data = data

    @classmethod
    def _load_lines(cls) -> list[str]:
        """Load the substructure SMARTS file, filtering out comments and blank lines.

        Returns:
            List of valid SMARTS lines from the file.

        Raises:
            FileNotFoundError: If the substructure file cannot be found.
            IOError: If there are issues reading the file.
        """
        logger.trace("Loading substructure data from site_substructures.smarts")

        try:
            with pkg_resources.open_text(
                "dimorphite_dl.smarts", "site_substructures.smarts"
            ) as f:
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

                logger.info("Loaded {} valid SMARTS patterns", valid_count)
                return lines

        except FileNotFoundError:
            logger.error("Could not find site_substructures.smarts file")
            raise
        except Exception as e:
            logger.error("Error reading substructure file: {}", str(e))
            raise IOError(f"Failed to read substructure file: {e}")

    @classmethod
    def _parse_substructure_line(cls, line: str) -> SubstructureDatum:
        """Parse a single line from the substructure data file.

        Args:
            line: Line from the substructure file

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
        return SubstructureDatum(name=name, smarts=smarts, pkas=pka_data, mol=mol)

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
    def _parse_pka_line(cls, line_parts: list[str]) -> list[PKaDatum]:
        if len(line_parts) % 3 != 0:
            logger.warning(
                "Invalid pKa data format, expected groups of 3, got {}", len(line_parts)
            )
            raise ValueError

        pka_data = []
        for i in range(0, len(line_parts), 3):
            try:
                idx_site = int(line_parts[i])
                mean = float(line_parts[i + 1])
                stdev = float(line_parts[i + 2])
                pka_data.append(PKaDatum(idx_site=idx_site, mean=mean, stdev=stdev))
            except (ValueError, IndexError) as e:
                logger.warning("Error parsing pKa data: {}", line_parts)
                raise ValueError from e
        return pka_data

    @classmethod
    def get_substructures(cls) -> Iterator[SubstructureDatum]:
        for substruct in cls._data:
            yield substruct
