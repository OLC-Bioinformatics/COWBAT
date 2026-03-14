# tests/test_fastqc.py

"""
Unit tests for FastQC functions in the Teacup COWBAT pipeline
"""

# Standard imports
import logging
import unittest
from unittest.mock import patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
from cowbat.fastqc import (
    fastqc,
    prepare_fastqc_call,
    get_reader_and_files,
    run_fastqc,
    _run_command,
    _log_fastqc_output,
    rename_fastqc_outputs,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class TestFastQC(unittest.TestCase):
    """
    Unit tests for FastQc analyses. Anything with queues is not covered, as
    the tests would just hang
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.log_file = 'tests/testdata/log.txt'
        cls.threads = 4
        cls.metadata = [
            CustomBox({
                'name': 'sample1',
                'general': {
                    'fastq_files': [
                        'sample1_R1.fastq.gz', 'sample1_R2.fastq.gz'
                    ],
                    'trimmed_fastq_files': [
                        'sample1_R1_trimmed.fastq.gz',
                        'sample1_R2_trimmed.fastq.gz'
                    ],
                    'trimmed_corrected_fastq_files': [],
                    'log_out': 'sample1_log_out.txt',
                    'log_err': 'sample1_log_err.txt',
                    'output_directory': 'tests/testdata/output'
                },
                'commands': {}
            })
        ]

    @patch('cowbat.fastqc.qc_queue')
    @patch('cowbat.fastqc.prepare_fastqc_call')
    def test_fastqc(self, mock_prepare_fastqc_call, mock_qc_queue):
        """
        Test that fastqc correctly prepares and processes FastQC calls.
        """
        mock_prepare_fastqc_call.return_value = (
            'fastqc_call', 'fastqc_reads'
        )
        mock_qc_queue.empty.return_value = True

        result = fastqc(
            self.metadata[0], 'trimmed', self.log_file, self.threads
        )

        mock_prepare_fastqc_call.assert_called_once()
        mock_qc_queue.put.assert_called_once()
        assert result == self.metadata[0]

    def test_prepare_fastqc_call(self):
        """
        Test that prepare_fastqc_call correctly prepares FastQC calls.
        """
        result = prepare_fastqc_call(
            self.metadata[0], 'trimmed', self.threads
        )

        assert result[0] != ''
        assert result[1] != ''

    def test_get_reader_and_files(self):
        """
        Test that get_reader_and_files correctly returns reader and files.
        """
        result = get_reader_and_files(self.metadata[0], 'fastq_files')

        assert result[0] == 'gunzip --to-stdout'
        assert result[1] == self.metadata[0].general.fastq_files

    @patch('cowbat.fastqc._run_command')
    @patch('cowbat.fastqc._log_fastqc_output')
    def test_run_fastqc(self, mock_log_fastqc_output, mock_run_command):
        """
        Test that run_fastqc correctly runs FastQC commands.
        """
        mock_run_command.return_value = ('stdout', 'stderr')

        result = run_fastqc(
            'system_call', 'fastqc_reads', self.log_file, 'trimmed',
            'tests/testdata/output', self.metadata[0]
        )

        mock_run_command.assert_called()
        mock_log_fastqc_output.assert_called()
        assert result == ('stdout stdout', 'stderr stderr')

    @patch('cowbat.fastqc.run_subprocess')
    def test_run_command(self, mock_run_subprocess):
        """
        Test that _run_command correctly runs a subprocess command.
        """
        mock_run_subprocess.return_value = ('stdout', 'stderr')

        result = _run_command('command')

        mock_run_subprocess.assert_called_once_with(command='command')
        assert result == ('stdout', 'stderr')

    @patch('cowbat.fastqc.write_to_log_file')
    def test_log_fastqc_output(self, mock_write_to_log_file):
        """
        Test that _log_fastqc_output correctly logs FastQC output.
        """
        _log_fastqc_output(
            'system_call', 'fastqc_reads', 'stdout', 'stderr',
            self.log_file, self.metadata[0]
        )

        assert mock_write_to_log_file.call_count == 3

    @patch('cowbat.fastqc.shutil.move')
    def test_rename_fastqc_outputs(self, mock_shutil_move):
        """
        Test that rename_fastqc_outputs correctly renames FastQC output files.
        """
        rename_fastqc_outputs(
            'trimmed', 'tests/testdata/output', self.metadata[0]
        )

        assert mock_shutil_move.call_count == 2
