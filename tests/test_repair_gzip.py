# tests/test_repair_gzip.py

"""
Unit tests for repairing gzip issues in the Teacup COWBAT pipeline
"""

# Standard imports
from unittest.mock import patch, MagicMock

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
from cowbat.repair_gzip import repair_gzip


class TestRepairGzip:
    """
    Unit tests for Teacup COWBAT repair gzip functions
    """
    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.log_file = 'tests/testdata/log.txt'
        cls.metadata = [
            CustomBox({
                'name': 'sample1',
                'general': {
                    'trimmed_corrected_fastq_files': [
                        'sample1_R1_trimmed_corrected.fastq.gz',
                        'sample1_R2_trimmed_corrected.fastq.gz'
                    ],
                    'output_directory': 'tests/testdata/output',
                    'log_out': 'sample1_log_out.txt',
                    'log_err': 'sample1_log_err.txt'
                }
            })
        ]

    @patch('cowbat.repair_gzip.run_subprocess')
    @patch('cowbat.repair_gzip.write_to_log_file')
    @patch('logging.info')
    @patch('logging.debug')
    def test_repair_gzip_success(
        self, mock_logging_debug, mock_logging_info,
        mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test that repair_gzip successfully processes samples.
        """
        print('Starting test_repair_gzip_success')
        # Mock run_subprocess to simulate successful execution
        mock_run_subprocess.return_value = ('stdout', 'stderr')

        # Call the function
        repair_gzip(self.log_file, self.metadata)

        # Assertions
        mock_logging_info.assert_any_call('Fixing BBMap gunzip issue')
        mock_logging_debug.assert_any_call('Processing sample: %s', 'sample1')
        mock_logging_debug.assert_any_call(
            'Running system call for forward read: %s',
            'gunzip -c sample1_R1_trimmed_corrected.fastq.gz > '
            'tests/testdata/output/sample1_tmp_R1.fastq && '
            'gzip -c tests/testdata/output/sample1_tmp_R1.fastq > '
            'sample1_R1_trimmed_corrected.fastq.gz && '
            'rm tests/testdata/output/sample1_tmp_R1.fastq'
        )
        mock_logging_debug.assert_any_call(
            'Running system call for reverse read: %s',
            'gunzip -c sample1_R2_trimmed_corrected.fastq.gz > '
            'tests/testdata/output/sample1_tmp_R2.fastq && '
            'gzip -c tests/testdata/output/sample1_tmp_R2.fastq > '
            'sample1_R2_trimmed_corrected.fastq.gz && '
            'rm tests/testdata/output/sample1_tmp_R2.fastq'
        )
        mock_write_to_log_file.assert_called()
        mock_run_subprocess.assert_called()
        print('Completed test_repair_gzip_success')

    @patch('cowbat.repair_gzip.run_subprocess')
    @patch('cowbat.repair_gzip.write_to_log_file')
    @patch('logging.warning')
    def test_repair_gzip_no_fastq_files(
        self, mock_logging_warning,
        mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test that repair_gzip handles missing FASTQ files.
        """
        print('Starting test_repair_gzip_no_fastq_files')
        # Modify metadata to have no FASTQ files
        self.metadata[0].general.trimmed_corrected_fastq_files = None

        # Call the function
        repair_gzip(self.log_file, self.metadata)

        # Assertions
        mock_logging_warning.assert_any_call(
            'No FASTQ files found for sample: %s', 'sample1'
        )
        mock_write_to_log_file.assert_not_called()
        mock_run_subprocess.assert_not_called()
        print('Completed test_repair_gzip_no_fastq_files')

    @patch('cowbat.repair_gzip.run_subprocess')
    @patch('cowbat.repair_gzip.write_to_log_file')
    @patch('logging.debug')
    def test_repair_gzip_attribute_error(
        self, mock_logging_debug, mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test that repair_gzip handles AttributeError.
        """
        print('Starting test_repair_gzip_attribute_error')
        # Modify metadata to raise AttributeError
        self.metadata[0].general = MagicMock()
        self.metadata[0].general.trimmed_corrected_fastq_files.side_effect = (
            AttributeError('attribute error')
        )

        # Call the function
        repair_gzip(self.log_file, self.metadata)

        # Assertions
        mock_logging_debug.assert_any_call('Processing sample: %s', 'sample1')
        mock_write_to_log_file.assert_not_called()
        mock_run_subprocess.assert_not_called()
        print('Completed test_repair_gzip_attribute_error')
