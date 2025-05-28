"""
Comprehensive unit tests for the SMILES processing library.

Tests cover all functionality including streaming, validation, error handling,
file processing, and edge cases for robust production use.
"""

import gzip
import os
import tempfile
from unittest.mock import patch

import pytest

from dimorphite_dl.io import (
    SMILESProcessor,
    SMILESRecord,
    SMILESStreamError,
    SMILESValidationError,
    process_smiles_file,
    stream_smiles,
)


class TestSMILESRecord:
    """Test SMILESRecord dataclass functionality."""

    def test_smiles_record_creation(self):
        """Test basic SMILESRecord creation."""
        record = SMILESRecord("CCO")
        assert record.smiles == "CCO"
        assert record.identifier == ""
        assert record.source_line is None
        assert record.metadata == {}

    def test_smiles_record_with_all_fields(self):
        """Test SMILESRecord with all fields populated."""
        metadata = {"property": "value"}
        record = SMILESRecord(
            smiles="CCO", identifier="ethanol", source_line=5, metadata=metadata
        )
        assert record.smiles == "CCO"
        assert record.identifier == "ethanol"
        assert record.source_line == 5
        assert record.metadata == metadata

    def test_smiles_record_metadata_initialization(self):
        """Test that metadata is properly initialized as empty dict."""
        record = SMILESRecord("CCO")
        assert record.metadata != dict()

        # Ensure it's a new dict each time
        record2 = SMILESRecord("CCC")
        assert record2.metadata != dict()

        record.metadata["test"] = "value"
        assert "test" not in record2.metadata


class TestSMILESProcessor:
    """Test SMILESProcessor functionality."""

    def test_processor_initialization_defaults(self):
        """Test processor initialization with default values."""
        processor = SMILESProcessor()
        assert processor.validate_smiles is True
        assert processor.skip_invalid is True
        assert processor.max_length == 10000
        assert processor.chunk_size == 1000
        assert processor.get_stats() == {"processed": 0, "skipped": 0, "errors": 0}

    def test_processor_initialization_custom(self):
        """Test processor initialization with custom values."""
        processor = SMILESProcessor(
            validate_smiles=False, skip_invalid=False, max_length=5000, chunk_size=500
        )
        assert processor.validate_smiles is False
        assert processor.skip_invalid is False
        assert processor.max_length == 5000
        assert processor.chunk_size == 500

    def test_single_smiles_string_input(self):
        """Test processing a single SMILES string."""
        processor = SMILESProcessor(validate_smiles=False)
        records = list(processor.stream("CCO"))

        assert len(records) == 1
        assert records[0].smiles == "CCO"
        assert records[0].source_line == 1
        assert processor.get_stats()["processed"] == 1

    def test_list_input(self):
        """Test processing a list of SMILES strings."""
        smiles_list = ["CCO", "CCC", "c1ccccc1"]
        processor = SMILESProcessor(validate_smiles=False)
        records = list(processor.stream(smiles_list))

        assert len(records) == 3
        assert [r.smiles for r in records] == smiles_list
        assert [r.source_line for r in records] == [1, 2, 3]
        assert processor.get_stats()["processed"] == 3

    def test_generator_input(self):
        """Test processing a generator of SMILES strings."""

        def smiles_generator():
            yield "CCO"
            yield "CCC"
            yield "c1ccccc1"

        processor = SMILESProcessor(validate_smiles=False)
        records = list(processor.stream(smiles_generator()))

        assert len(records) == 3
        assert processor.get_stats()["processed"] == 3

    def test_empty_input_handling(self):
        """Test handling of empty inputs."""
        processor = SMILESProcessor(validate_smiles=False)

        # Empty list
        records = list(processor.stream([]))
        assert len(records) == 0

        # List with empty strings
        records = list(processor.stream(["", "  ", "CCO"]))
        assert len(records) == 1
        assert records[0].smiles == "CCO"

    def test_whitespace_stripping(self):
        """Test that whitespace is properly stripped."""
        processor = SMILESProcessor(validate_smiles=False)
        records = list(processor.stream(["  CCO  ", "\tCCC\n"]))

        assert len(records) == 2
        assert records[0].smiles == "CCO"
        assert records[1].smiles == "CCC"

    def test_max_length_validation(self):
        """Test maximum length validation."""
        processor = SMILESProcessor(max_length=5, validate_smiles=False)
        records = list(processor.stream(["CCO", "CCCCCCCCCC"]))

        assert len(records) == 1
        assert records[0].smiles == "CCO"
        stats = processor.get_stats()
        assert stats["processed"] == 1
        assert stats["skipped"] == 1
        assert stats["errors"] == 1

    def test_max_length_none(self):
        """Test that max_length=None allows any length."""
        very_long_smiles = "C" * 20000
        processor = SMILESProcessor(max_length=None, validate_smiles=False)
        records = list(processor.stream([very_long_smiles]))

        assert len(records) == 1
        assert records[0].smiles == very_long_smiles

    @patch("dimorphite_dl.io.rdMolStandardize.ValidateSmiles")
    def test_smiles_validation_valid(self, mock_validate):
        """Test SMILES validation with valid molecules."""
        mock_validate.return_value = None  # No exception = valid

        processor = SMILESProcessor(validate_smiles=True)
        records = list(processor.stream(["CCO"]))

        assert len(records) == 1
        assert records[0].smiles == "CCO"
        mock_validate.assert_called_once_with("CCO")

    @patch("dimorphite_dl.io.rdMolStandardize.ValidateSmiles")
    def test_smiles_validation_invalid_skip(self, mock_validate):
        """Test SMILES validation with invalid molecules (skip mode)."""
        mock_validate.side_effect = Exception("Invalid SMILES")

        processor = SMILESProcessor(validate_smiles=True, skip_invalid=True)
        records = list(processor.stream(["invalid_smiles"]))

        assert len(records) == 0
        stats = processor.get_stats()
        assert stats["skipped"] == 1
        assert stats["errors"] == 1

    @patch("dimorphite_dl.io.rdMolStandardize.ValidateSmiles")
    def test_smiles_validation_invalid_fail(self, mock_validate):
        """Test SMILES validation with invalid molecules (fail mode)."""
        mock_validate.side_effect = Exception("Invalid SMILES")

        processor = SMILESProcessor(validate_smiles=True, skip_invalid=False)

        with pytest.raises(SMILESValidationError):
            list(processor.stream(["invalid_smiles"]))

    def test_non_string_items_in_iterable(self):
        """Test handling of non-string items in iterable."""
        processor = SMILESProcessor(validate_smiles=False)
        mixed_input = ["CCO", 123, "CCC", None]

        records = list(processor.stream(mixed_input))  # type: ignore
        assert len(records) == 2
        assert [r.smiles for r in records] == ["CCO", "CCC"]

        stats = processor.get_stats()
        assert stats["errors"] == 2

    def test_unsupported_input_type(self):
        """Test handling of unsupported input types."""
        processor = SMILESProcessor()

        list(processor.stream(123))  # type: ignore

    def test_stats_reset_between_calls(self):
        """Test that stats are reset between stream calls."""
        processor = SMILESProcessor(validate_smiles=False)

        list(processor.stream(["CCO", "CCC"]))
        stats1 = processor.get_stats()

        list(processor.stream(["c1ccccc1"]))
        stats2 = processor.get_stats()

        assert stats1["processed"] == 2
        assert stats2["processed"] == 1  # Reset for new call

    def test_batch_streaming(self):
        """Test batch streaming functionality."""
        smiles_list = ["CCO", "CCC", "c1ccccc1", "CCCC", "CCCCC"]
        processor = SMILESProcessor(validate_smiles=False)

        batches = list(processor.stream_batches(smiles_list, batch_size=2))

        assert len(batches) == 3
        assert len(batches[0]) == 2
        assert len(batches[1]) == 2
        assert len(batches[2]) == 1  # Remainder

        # Check SMILES are correct
        all_smiles = []
        for batch in batches:
            all_smiles.extend([r.smiles for r in batch])
        assert all_smiles == smiles_list

    def test_batch_streaming_default_size(self):
        """Test batch streaming with default chunk size."""
        processor = SMILESProcessor(chunk_size=3, validate_smiles=False)
        smiles_list = ["CCO"] * 7

        batches = list(processor.stream_batches(smiles_list))

        assert len(batches) == 3
        assert len(batches[0]) == 3
        assert len(batches[1]) == 3
        assert len(batches[2]) == 1


class TestFileProcessing:
    """Test file processing functionality."""

    def test_is_file_path_detection(self):
        """Test file path detection logic."""
        processor = SMILESProcessor()

        # Clear file paths
        assert processor._is_file_path("/path/to/file.smi")
        assert processor._is_file_path("file.csv")
        assert processor._is_file_path("data.txt")
        assert processor._is_file_path("molecules.smiles")

        # Clear SMILES strings
        assert not processor._is_file_path("CCO")
        assert not processor._is_file_path("c1ccccc1")

        # Very long strings should not be file paths
        long_string = "C" * 2000
        assert not processor._is_file_path(long_string)

    def test_file_not_found(self):
        """Test handling of non-existent files."""
        import loguru
        from loguru import logger

        processor = SMILESProcessor()  # skip_invalid=True by default

        # Capture loguru logs
        import io

        log_stream = io.StringIO()

        # Add a sink to capture logs
        sink_id = logger.add(log_stream, format="{message}")

        try:
            records = list(processor.stream("nonexistent_file.smi"))

            # Should return empty results but log the error
            assert len(records) == 0

            # Check that error was logged
            log_output = log_stream.getvalue()
            assert "File not found: nonexistent_file.smi" in log_output

        finally:
            # Clean up the logger sink
            logger.remove(sink_id)

    def test_text_file_processing(self):
        """Test processing of plain text SMILES files."""
        content = "CCO\nCCC ethane\nc1ccccc1 benzene\n# comment\n\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            processor = SMILESProcessor(validate_smiles=False)
            records = list(processor.stream(temp_path))

            assert len(records) == 3
            assert records[0].smiles == "CCO"
            assert records[0].identifier == ""
            assert records[1].smiles == "CCC"
            assert records[1].identifier == "ethane"
            assert records[2].smiles == "c1ccccc1"
            assert records[2].identifier == "benzene"

            # Check line numbers
            assert records[0].source_line == 1
            assert records[1].source_line == 2
            assert records[2].source_line == 3

        finally:
            os.unlink(temp_path)

    def test_gzipped_file_processing(self):
        """Test processing of gzipped files."""
        content = "CCO\nCCC\n"

        with tempfile.NamedTemporaryFile(suffix=".smi.gz", delete=False) as f:
            temp_path = f.name

        try:
            with gzip.open(temp_path, "wt") as f:
                f.write(content)

            processor = SMILESProcessor(validate_smiles=False)
            records = list(processor.stream(temp_path))

            assert len(records) == 2
            assert records[0].smiles == "CCO"
            assert records[1].smiles == "CCC"

        finally:
            os.unlink(temp_path)

    def test_unknown_file_format(self):
        """Test handling of unknown file formats."""
        content = "CCO\nCCC\n"

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".unknown", delete=False
        ) as f:
            f.write(content)
            temp_path = f.name

        try:
            processor = SMILESProcessor(validate_smiles=False)
            records = list(processor.stream(temp_path))

            # Should treat as text file
            assert len(records) == 2
            assert records[0].smiles == "CCO"
            assert records[1].smiles == "CCC"

        finally:
            os.unlink(temp_path)

    def test_file_encoding_error_handling(self):
        """Test handling of files with encoding issues."""
        # Create file with binary data that will cause encoding issues
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".smi", delete=False) as f:
            f.write(b"CCO\n\xff\xfe\nCCC\n")  # Invalid UTF-8 bytes
            temp_path = f.name

        try:
            processor = SMILESProcessor(validate_smiles=False)
            records = list(processor.stream(temp_path))

            # Should handle encoding errors gracefully
            assert len(records) >= 2  # At least CCO and CCC should be processed

        finally:
            os.unlink(temp_path)

    @patch("dimorphite_dl.io.open")
    def test_file_read_error(self, mock_open_func):
        """Test handling of file read errors."""
        mock_open_func.side_effect = PermissionError("Permission denied")

        processor = SMILESProcessor(skip_invalid=False)

        with pytest.raises(SMILESStreamError):
            list(processor.stream("test.smi"))


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_stream_smiles_function(self):
        """Test stream_smiles convenience function."""
        smiles_list = ["CCO", "CCC"]
        records = list(stream_smiles(smiles_list, validate_smiles=False))

        assert len(records) == 2
        assert records[0].smiles == "CCO"
        assert records[1].smiles == "CCC"

    def test_process_smiles_file_function(self):
        """Test process_smiles_file convenience function."""
        content = "CCO\nCCC\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            records = list(process_smiles_file(temp_path, validate_smiles=False))

            assert len(records) == 2
            assert records[0].smiles == "CCO"
            assert records[1].smiles == "CCC"

        finally:
            os.unlink(temp_path)


class TestErrorHandling:
    """Test comprehensive error handling scenarios."""

    def test_skip_invalid_true(self):
        """Test skip_invalid=True behavior."""
        processor = SMILESProcessor(
            validate_smiles=False, skip_invalid=True, max_length=3
        )

        # Mix of valid and invalid (too long) SMILES
        smiles_list = ["CCO", "CCCCCCCCCC", "CCC"]
        records = list(processor.stream(smiles_list))

        assert len(records) == 2
        assert records[0].smiles == "CCO"
        assert records[1].smiles == "CCC"

        stats = processor.get_stats()
        assert stats["processed"] == 2
        assert stats["skipped"] == 1
        assert stats["errors"] == 1

    def test_skip_invalid_false(self):
        """Test skip_invalid=False behavior."""
        processor = SMILESProcessor(
            validate_smiles=False, skip_invalid=False, max_length=3
        )

        with pytest.raises(SMILESValidationError):
            list(processor.stream(["CCO", "CCCCCCCCCC"]))

    def test_exception_during_record_creation(self):
        """Test handling of unexpected exceptions during record creation."""
        processor = SMILESProcessor(validate_smiles=False)

        # Mock _validate_smiles_syntax to raise unexpected exception
        with patch.object(
            processor,
            "_validate_smiles_syntax",
            side_effect=RuntimeError("Unexpected error"),
        ):
            processor.validate_smiles = True
            records = list(processor.stream(["CCO"]))

            assert len(records) == 0
            stats = processor.get_stats()
            assert stats["errors"] == 1


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_smiles_handling(self):
        """Test handling of empty and whitespace-only SMILES."""
        processor = SMILESProcessor(validate_smiles=False)
        inputs = ["", "  ", "\t", "\n", "CCO", "   "]

        records = list(processor.stream(inputs))

        assert len(records) == 1
        assert records[0].smiles == "CCO"

    def test_comment_lines_in_file(self):
        """Test that comment lines are properly skipped."""
        content = "# This is a comment\nCCO\n# Another comment\nCCC\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            processor = SMILESProcessor(validate_smiles=False)
            records = list(processor.stream(temp_path))

            assert len(records) == 2
            assert records[0].smiles == "CCO"
            assert records[1].smiles == "CCC"

        finally:
            os.unlink(temp_path)

    def test_multiline_input_with_spaces(self):
        """Test handling of complex whitespace scenarios."""
        inputs = ["CCO  ethanol", " CCC\tpropane ", "c1ccccc1"]
        processor = SMILESProcessor(validate_smiles=False)

        records = list(processor.stream(inputs))

        assert len(records) == 3
        assert records[0].smiles == "CCO"
        assert records[0].identifier == "ethanol"
        assert records[1].smiles == "CCC"
        assert records[1].identifier == "propane"
        assert records[2].smiles == "c1ccccc1"
        assert records[2].identifier == ""

    def test_very_large_batch_size(self):
        """Test batch processing with very large batch size."""
        smiles_list = ["CCO"] * 5
        processor = SMILESProcessor(validate_smiles=False)

        batches = list(processor.stream_batches(smiles_list, batch_size=1000))

        assert len(batches) == 1
        assert len(batches[0]) == 5

    def test_batch_size_one(self):
        """Test batch processing with batch size of 1."""
        smiles_list = ["CCO", "CCC", "c1ccccc1"]
        processor = SMILESProcessor(validate_smiles=False)

        batches = list(processor.stream_batches(smiles_list, batch_size=1))

        assert len(batches) == 3
        for batch in batches:
            assert len(batch) == 1

    def test_metadata_preservation(self):
        """Test that metadata is properly preserved and isolated."""
        processor = SMILESProcessor(validate_smiles=False)

        record1 = processor._create_record("CCO", metadata={"prop": "value1"})
        record2 = processor._create_record("CCC", metadata={"prop": "value2"})

        assert record1.metadata["prop"] == "value1"  # type: ignore
        assert record2.metadata["prop"] == "value2"  # type: ignore

        # Modifying one shouldn't affect the other
        record1.metadata["new"] = "test"  # type: ignore
        assert "new" not in record2.metadata  # type: ignore


class TestIntegration:
    """Integration tests combining multiple features."""

    def test_large_file_simulation(self):
        """Test processing a simulated large file."""
        # Create a larger test file
        content = "\n".join([f"{'C' * (i % 10 + 1)} mol_{i}" for i in range(1000)])

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            processor = SMILESProcessor(validate_smiles=False)

            # Process in batches
            total_processed = 0
            for batch in processor.stream_batches(temp_path, batch_size=100):
                total_processed += len(batch)
                assert len(batch) <= 100

            assert total_processed == 1000

        finally:
            os.unlink(temp_path)

    def test_mixed_valid_invalid_with_validation(self):
        """Test processing mixed valid/invalid SMILES with validation enabled."""
        # Note: This test would need actual RDKit validation,
        # so we'll mock it for testing purposes

        with patch("dimorphite_dl.io.rdMolStandardize.ValidateSmiles") as mock_validate:
            # First two are valid, third is invalid
            mock_validate.side_effect = [None, None, Exception("Invalid")]

            processor = SMILESProcessor(validate_smiles=True, skip_invalid=True)
            smiles_list = ["CCO", "CCC", "invalid"]

            records = list(processor.stream(smiles_list))

            assert len(records) == 2
            assert records[0].smiles == "CCO"
            assert records[1].smiles == "CCC"

            stats = processor.get_stats()
            assert stats["processed"] == 2
            assert stats["skipped"] == 1

    def test_end_to_end_file_processing(self):
        """Test complete end-to-end file processing workflow."""
        # Create test file with various scenarios
        content = """# Test SMILES file
CCO ethanol
CCC propane
# Another comment

CCCC butane
invalid_smiles_here should_be_skipped
c1ccccc1 benzene
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".smi", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            processor = SMILESProcessor(
                validate_smiles=False,  # Skip RDKit validation for test
                skip_invalid=True,
                max_length=100,
            )

            records = list(processor.stream(temp_path))

            expected_smiles = ["CCO", "CCC", "CCCC", "invalid_smiles_here", "c1ccccc1"]
            expected_ids = [
                "ethanol",
                "propane",
                "butane",
                "should_be_skipped",
                "benzene",
            ]

            assert len(records) == 5
            for i, record in enumerate(records):
                assert record.smiles == expected_smiles[i]
                assert record.identifier == expected_ids[i]
                assert record.source_line is not None
                assert isinstance(record.metadata, dict)

            stats = processor.get_stats()
            assert stats["processed"] == 5
            assert stats["errors"] == 0

        finally:
            os.unlink(temp_path)
