"""
Robust, memory-efficient SMILES string handling library.

Provides unified streaming interface for processing SMILES from various sources
with comprehensive error handling, validation, and memory optimization.
"""

from typing import Any, TextIO

import gzip
import os
import pathlib
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field

from loguru import logger
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass
class SMILESRecord:
    """Container for a SMILES string with metadata."""

    smiles: str
    identifier: str = ""
    source_line: int | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


class SMILESValidationError(Exception):
    """Raised when SMILES validation fails."""

    pass


class SMILESStreamError(Exception):
    """Raised when streaming encounters an error."""

    pass


class SMILESProcessor:
    """
    Memory-efficient SMILES string processor with robust error handling.

    Handles various input formats and provides streaming interface for
    processing large datasets without memory overflow.
    """

    def __init__(
        self,
        validate_smiles: bool = True,
        skip_invalid: bool = True,
        max_length: int | None = 10000,
        chunk_size: int = 1000,
    ):
        """
        Initialize SMILES processor.

        Args:
            validate_smiles: Whether to validate SMILES syntax
            skip_invalid: Skip invalid SMILES instead of raising errors
            max_length: Maximum allowed SMILES length (None for no limit)
            chunk_size: Batch size for processing operations
        """
        self.validate_smiles = validate_smiles
        self.skip_invalid = skip_invalid
        self.max_length = max_length
        self.chunk_size = chunk_size
        self._stats: dict[str, int] = {"processed": 0, "skipped": 0, "errors": 0}

    def stream(
        self, input_data: str | Iterable[str] | Iterator[str]
    ) -> Iterator[SMILESRecord]:
        """
        Stream SMILES records from various input types.

        Args:
            input_data: File path, single SMILES, or iterable of SMILES

        Yields:
            SMILESRecord: Validated SMILES records with metadata

        Raises:
            SMILESStreamError: If input cannot be processed
        """
        self._reset_stats()

        try:
            if isinstance(input_data, str):
                yield from self._handle_string_input(input_data)
            elif hasattr(input_data, "__iter__"):
                yield from self._handle_iterable_input(input_data)
            else:
                raise SMILESStreamError(f"Unsupported input type: {type(input_data)}")

        except Exception as e:
            logger.error(f"Error streaming SMILES: {e}")
            if not self.skip_invalid:
                if isinstance(e, SMILESValidationError):
                    raise e
                else:
                    raise SMILESStreamError(f"Failed to process input: {e}") from e

    def stream_batches(
        self, input_data: str | Iterable[str], batch_size: int | None = None
    ) -> Iterator[list[SMILESRecord]]:
        """
        Stream SMILES records in batches for efficient processing.

        Args:
            input_data: Input source
            batch_size: Size of each batch (uses instance default if None)

        Yields:
            Batches of SMILES records
        """
        batch_size = batch_size or self.chunk_size
        batch = []

        for record in self.stream(input_data):
            batch.append(record)
            if len(batch) >= batch_size:
                yield batch
                batch = []

        if batch:  # Yield remaining records
            yield batch

    def _handle_string_input(self, input_str: str) -> Iterator[SMILESRecord]:
        """Handle string input - either file path or single SMILES."""
        if self._is_file_path(input_str):
            yield from self._stream_from_file(input_str)
        else:
            # Single SMILES string
            record = self._create_record(input_str, source_line=1)
            if record:
                yield record

    def _handle_iterable_input(self, iterable: Iterable[str]) -> Iterator[SMILESRecord]:
        """Handle iterable input (list, generator, etc.).

        This will skip empty lines.
        """
        logger.debug("Handling iterable input of {}", iterable)
        for line_num, line in enumerate(iterable, 1):
            if isinstance(line, str):
                line = line.strip()
                line_split = line.split()
                if len(line_split) > 2:
                    logger.warning(
                        f"Lines can only contain a smiles string and identifier, but we were given {line}"
                    )
                    raise ValueError(
                        "Line contains more than two items (smiles and identifier)"
                    )
                if len(line_split) == 0:
                    continue
                smiles = line_split[0]
                if len(line_split) == 2:
                    identifier = line_split[1]
                else:
                    identifier = ""
                record = self._create_record(smiles, identifier, source_line=line_num)
                if record:
                    yield record
            else:
                self._handle_error(
                    f"Non-string item at position {line_num}: {type(line)}"
                )

    def _stream_from_file(self, filepath: str) -> Iterator[SMILESRecord]:
        """Stream SMILES from file with format auto-detection."""
        logger.debug("Streaming from {}", filepath)
        path = pathlib.Path(filepath)

        if not path.exists():
            raise SMILESStreamError(f"File not found: {filepath}")

        # Handle compressed files
        open_func = gzip.open if path.suffix == ".gz" else open
        mode = "rt" if path.suffix == ".gz" else "r"

        try:
            with open_func(filepath, mode, encoding="utf-8", errors="replace") as f:
                yield from self._stream_from_file_object(f, path)
        except Exception as e:
            raise SMILESStreamError(f"Error reading file {filepath}: {e}") from e

    def _stream_from_file_object(
        self, file_obj: TextIO, path: pathlib.Path
    ) -> Iterator[SMILESRecord]:
        """Stream from file object based on file extension."""
        suffix = path.suffix.lower().replace(".gz", "")

        if suffix in {".smiles", ".smi", ".txt", ""}:
            yield from self._stream_from_text(file_obj)
        else:
            logger.warning(f"Unknown file format {suffix}, treating as text")
            yield from self._stream_from_text(file_obj)

    def _stream_from_text(self, file_obj: TextIO) -> Iterator[SMILESRecord]:
        """Stream SMILES from plain text file."""
        for line_num, line in enumerate(file_obj, 1):
            line = line.strip()
            if line and not line.startswith("#"):
                # Handle multi-column format (SMILES ID)
                parts = line.split()
                smiles = parts[0]
                identifier = parts[1] if len(parts) > 1 else ""

                record = self._create_record(
                    smiles, identifier=identifier, source_line=line_num
                )
                if record:
                    yield record

    def _create_record(
        self,
        smiles: str,
        identifier: str = "",
        source_line: int | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> SMILESRecord | None:
        """Create and validate a SMILES record."""
        smiles = smiles.strip()

        if not smiles:
            return None

        try:
            # Length validation
            if self.max_length and len(smiles) > self.max_length:
                self._handle_error(
                    f"SMILES too long ({len(smiles)} > {self.max_length}): {smiles[:50]}..."
                )
                return None

            # Basic syntax validation
            if self.validate_smiles:
                if not self._validate_smiles_syntax(smiles):
                    self._handle_error(f"Invalid SMILES syntax: {smiles}")
                    return None

            self._stats["processed"] += 1
            return SMILESRecord(
                smiles=smiles,
                identifier=identifier,
                source_line=source_line,
                metadata=metadata or {},
            )

        except Exception as e:
            self._handle_error(f"Error creating record for '{smiles}': {e}")
            return None

    def _validate_smiles_syntax(self, smiles: str) -> bool:
        """SMILES syntax validation using RDKit."""
        logger.info("Processing {}", smiles)
        try:
            rdMolStandardize.ValidateSmiles(smiles)
            logger.debug("SMILES is valid")
            return True
        except Exception:
            logger.info("SMILES is NOT valid")
            return False

    def _is_file_path(self, s: str) -> bool:
        """Check if string is likely a file path."""
        # Don't treat very long strings as file paths
        if len(s) > 1000:
            return False

        return bool(
            os.path.exists(s)
            or os.path.sep in s
            or (os.path.altsep and os.path.altsep in s)
            or s.endswith((".smiles", ".smi", ".txt", ".csv", ".sdf", ".gz"))
        )

    def _handle_error(self, message: str):
        """Handle errors based on skip_invalid setting."""
        self._stats["errors"] += 1
        if self.skip_invalid:
            self._stats["skipped"] += 1
            logger.warning(message)
        else:
            raise SMILESValidationError(message)

    def _reset_stats(self):
        """Reset processing statistics."""
        self._stats = {"processed": 0, "skipped": 0, "errors": 0}

    def get_stats(self) -> dict[str, int]:
        """Get processing statistics."""
        return self._stats.copy()


# Convenience functions
def stream_smiles(input_data: str | Iterable[str], **kwargs) -> Iterator[SMILESRecord]:
    """Convenience function for streaming SMILES."""
    processor = SMILESProcessor(**kwargs)
    yield from processor.stream(input_data)


def process_smiles_file(filepath: str, **kwargs) -> Iterator[SMILESRecord]:
    """Convenience function for processing SMILES files."""
    processor = SMILESProcessor(**kwargs)
    yield from processor.stream(filepath)
