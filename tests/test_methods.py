#!/usr/env/python3

"""
Unit tests for Teacup COWBAT
"""

# Standard imports
import logging
import os
import subprocess
from unittest.mock import MagicMock, patch

# Third-party imports
from cowbat.methods import (
    check_programs,
    initialize_logging,
    read_checkpoint,
    sample_metadata,
    tilde_expand,
    write_checkpoint
)
from cowbat.teacup_version import __version__


class TestTeacup:
    """
    Unit tests for the Teacup COWBAT pipeline
    """
    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize the metadata before any tests run.
        """
        cls.metadata = []
        cls.checkpoint_file = 'tests/testdata/.checkpoint'
        cls.nonexistent_file = 'tests/testdata/nonexistent_checkpoint'
        cls.home_dir = os.path.expanduser('~')
        cls.test_dir = 'tests/testdata'
        cls.log_file = 'tests/testdata/test_log.log'
        cls.error_log_file = 'tests/testdata/test_error_log.log'
        cls.logging_level = 'DEBUG'
        cls.sample_sheet_path = os.path.join(cls.test_dir, 'SampleSheet.csv')
        cls.error_logger = MagicMock()
        cls.programs = ['python', 'nonexistent_program']

    def test_read_checkpoint_file_exists(self):
        """
        Test that read_checkpoint correctly reads the last successful step
        from an existing checkpoint file.
        """
        # Create a checkpoint file with a sample step
        with open(self.checkpoint_file, 'w', encoding='utf-8') as f:
            f.write('step1')

        # Call the function and check the result
        result = read_checkpoint(self.checkpoint_file)
        assert result == 'step1'

        # Clean up the checkpoint file
        os.remove(self.checkpoint_file)

    def test_read_checkpoint_file_does_not_exist(self):
        """
        Test that read_checkpoint returns an empty string when the checkpoint
        file does not exist.
        """
        # Ensure the nonexistent file does not exist
        if os.path.exists(self.nonexistent_file):
            os.remove(self.nonexistent_file)

        # Call the function and check the result
        result = read_checkpoint(self.nonexistent_file)
        assert result == ''

    def test_write_checkpoint(self):
        """
        Test that write_checkpoint correctly writes the current step to the
        checkpoint file.
        """
        step = 'step1'

        # Call the function to write the step to the checkpoint file
        write_checkpoint(self.checkpoint_file, step)

        # Verify that the step was written correctly
        with open(self.checkpoint_file, 'r', encoding='utf-8') as f:
            content = f.read().strip()

        assert content == step

        # Clean up the checkpoint file
        os.remove(self.checkpoint_file)

    def test_tilde_expand_with_tilde(self):
        """
        Test that tilde_expand correctly expands the tilde to the full user
        directory path.
        """
        path_with_tilde = '~/my_dir'
        expected_path = os.path.join(self.home_dir, 'my_dir')

        # Call the function and check the result
        result = tilde_expand(path_with_tilde)
        assert result == expected_path

    def test_tilde_expand_without_tilde(self):
        """
        Test that tilde_expand returns the absolute path when no tilde is
        present.
        """
        relative_path = 'tests/testdata'
        expected_path = os.path.abspath(relative_path)

        # Call the function and check the result
        result = tilde_expand(relative_path)
        assert result == expected_path

    def test_initialize_logging_creates_log_files(self):
        """
        Test that initialize_logging creates the log files and sets up the
        loggers correctly.
        """
        # Call the function to initialize logging
        logger, error_logger = initialize_logging(
            self.error_log_file, self.log_file, self.logging_level
        )

        # Check that the log files are created
        assert os.path.exists(self.log_file)
        assert os.path.exists(self.error_log_file)

        # Check that the loggers are set up correctly
        assert logger.level == logging.DEBUG
        assert error_logger.level == logging.ERROR

        # Clean up the log files
        os.remove(self.log_file)
        os.remove(self.error_log_file)

    def test_initialize_logging_writes_to_log_files(self):
        """
        Test that initialize_logging writes log messages to the log files.
        """
        # Call the function to initialize logging
        logger, error_logger = initialize_logging(
            self.error_log_file, self.log_file, self.logging_level
        )

        # Write a log message
        logger.debug('This is a debug message')
        error_logger.error('This is an error message')

        # Check that the log messages are written to the log files
        with open(self.log_file, 'r', encoding='utf-8') as f:
            log_content = f.read()
            assert 'This is a debug message' in log_content

        with open(self.error_log_file, 'r', encoding='utf-8') as f:
            error_log_content = f.read()
            assert 'This is an error message' in error_log_content

        # Clean up the log files
        os.remove(self.log_file)
        os.remove(self.error_log_file)

    @patch('cowbat.methods.determine_and_parse_sample_sheet')
    @patch('cowbat.methods.process_sample')
    @patch('cowbat.methods.FastqMover')
    @patch('cowbat.methods.write_metadata_to_file')
    def test_sample_metadata_with_sample_sheet(
        self, mock_write_metadata_to_file, mock_fastq_mover,
        mock_process_sample, mock_determine_and_parse_sample_sheet
    ):
        """
        Test that sample_metadata processes the sample sheet correctly when
        it is present.
        """
        # Create a mock sample sheet
        with open(self.sample_sheet_path, 'w', encoding='utf-8') as f:
            f.write('SampleID,SampleName\n1,Sample1\n2,Sample2')

        # Mock the return values of the patched functions
        mock_determine_and_parse_sample_sheet.return_value = {
            'SampleID': [1, 2], 'SampleName': ['Sample1', 'Sample2']
        }
        mock_process_sample.return_value = [
            {'id': 1, 'name': 'Sample1'}, {'id': 2, 'name': 'Sample2'}
        ]

        # Call the function
        metadata = sample_metadata(self.error_logger, [], self.test_dir)

        # Assertions
        mock_determine_and_parse_sample_sheet.assert_called_once_with(
            file_path=self.sample_sheet_path
        )
        mock_process_sample.assert_called_once()
        mock_fastq_mover.assert_called_once_with(
            metadata=metadata, path=self.test_dir
        )
        mock_write_metadata_to_file.assert_called_once_with(
            error_logger=self.error_logger, metadata=metadata
        )

        # Clean up the sample sheet
        os.remove(self.sample_sheet_path)

    @patch('cowbat.methods.basic')
    @patch('cowbat.methods.FastqMover')
    @patch('cowbat.methods.write_metadata_to_file')
    def test_sample_metadata_without_sample_sheet(
        self, mock_write_metadata_to_file, mock_fastq_mover, mock_basic
    ):
        """
        Test that sample_metadata processes correctly when the sample sheet
        is not present.
        """
        # Ensure the sample sheet does not exist
        if os.path.exists(self.sample_sheet_path):
            os.remove(self.sample_sheet_path)

        # Mock the return values of the patched functions
        mock_basic.return_value = [
            {'id': 1, 'name': 'Sample1'}, {'id': 2, 'name': 'Sample2'}
        ]

        # Call the function
        metadata = sample_metadata(self.error_logger, [], self.test_dir)

        # Assertions
        mock_basic.assert_called_once_with(
            commit=__version__, sequence_path=self.test_dir
        )
        mock_fastq_mover.assert_called_once_with(
            metadata=metadata, path=self.test_dir
        )
        mock_write_metadata_to_file.assert_called_once_with(
            error_logger=self.error_logger, metadata=metadata
        )

    @patch('cowbat.methods.subprocess.run')
    @patch('cowbat.methods.logging')
    def test_check_programs_all_installed(
        self, mock_logging, mock_subprocess_run
    ):
        """
        Test that check_programs correctly identifies all programs as
        installed.
        """
        # Mock subprocess.run to simulate installed programs
        mock_subprocess_run.return_value = MagicMock(returncode=0)

        # Call the function
        result = check_programs(['python'])

        # Assertions
        mock_subprocess_run.assert_called_once_with(
            ['python', '--version'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
        mock_logging.info.assert_any_call("%s is installed.", 'python')
        mock_logging.info.assert_any_call(
            "All required programs are installed."
        )
        assert result is None

    @patch('cowbat.methods.subprocess.run')
    @patch('cowbat.methods.logging')
    def test_check_programs_some_missing(
        self, mock_logging, mock_subprocess_run
    ):
        """
        Test that check_programs correctly identifies some programs as missing.
        """
        # Mock subprocess.run to simulate one installed and one missing program
        def mock_run_side_effect(cmd, **_):
            if cmd[0] == 'python':
                return MagicMock(returncode=0)
            raise FileNotFoundError

        mock_subprocess_run.side_effect = mock_run_side_effect

        # Call the function
        result = check_programs(self.programs)

        # Assertions
        mock_logging.info.assert_any_call("%s is installed.", 'python')
        mock_logging.warning.assert_any_call(
            "%s is missing.", 'nonexistent_program'
        )
        mock_logging.error.assert_any_call(
            "The following programs are missing:"
        )
        mock_logging.error.assert_any_call('nonexistent_program')
        assert result == ['nonexistent_program']

    @patch('cowbat.methods.logging')
    def test_check_programs_no_programs_provided(self, mock_logging):
        """
        Test that check_programs logs an error when no programs are provided.
        """
        # Call the function
        result = check_programs([])

        # Assertions
        mock_logging.error.assert_any_call("No programs provided to check.")
        assert result is None
