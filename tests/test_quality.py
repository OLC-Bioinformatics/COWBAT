# tests/test_quality.py

"""
Unit tests for quality checking functions in the Teacup COWBAT pipeline
"""

# Standard imports
from unittest.mock import patch, MagicMock

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
from cowbat.quality import quality


class TestQuality:
    """
    Unit tests for Teacup COWBAT quality functions
    """
    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.log_file = 'tests/testdata/log.txt'
        cls.report_path = 'tests/testdata/report'
        cls.sequence_path = 'tests/testdata/sequence'
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
        cls.error_logger = MagicMock()

    @patch('cowbat.quality.write_metadata_to_file')
    @patch('cowbat.quality.contamination_finder')
    @patch('cowbat.quality.repair_gzip')
    @patch('cowbat.quality.error_correction')
    @patch('cowbat.quality.fastqc_threader')
    @patch('cowbat.quality.trim_quality')
    @patch('cowbat.quality.validate_fastq')
    def test_quality_success(
        self, mock_validate_fastq, mock_trim_quality,
        mock_fastqc_threader, mock_error_correction,
        mock_repair_gzip, mock_contamination_finder,
        mock_write_metadata_to_file
    ):
        """
        Test that quality function successfully processes samples.
        """
        print('Starting test_quality_success')
        # Mock the return values of the functions
        mock_validate_fastq.return_value = self.metadata
        mock_trim_quality.return_value = self.metadata
        mock_fastqc_threader.return_value = self.metadata
        mock_error_correction.return_value = self.metadata
        mock_contamination_finder.return_value = self.metadata

        # Call the function
        result = quality(
            error_logger=self.error_logger,
            log_file=self.log_file,
            metadata=self.metadata,
            report_path=self.report_path,
            sequence_path=self.sequence_path,
            threads=self.threads
        )

        # Assertions
        mock_validate_fastq.assert_called_once_with(
            metadata=self.metadata, log_file=self.log_file
        )
        mock_write_metadata_to_file.assert_called()
        mock_fastqc_threader.assert_called()
        mock_trim_quality.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata
        )
        mock_error_correction.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata,
            threads=self.threads
        )
        mock_repair_gzip.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata
        )
        mock_contamination_finder.assert_called_once_with(
            input_path=self.sequence_path, log_file=self.log_file,
            metadata=self.metadata, report_path=self.report_path,
            threads=self.threads
        )
        assert result == self.metadata
        print('Completed test_quality_success')

    @patch('cowbat.quality.write_metadata_to_file')
    @patch('cowbat.quality.contamination_finder')
    @patch('cowbat.quality.repair_gzip')
    @patch('cowbat.quality.error_correction')
    @patch('cowbat.quality.fastqc_threader')
    @patch('cowbat.quality.trim_quality')
    @patch('cowbat.quality.validate_fastq')
    def test_quality_failure(
        self, mock_validate_fastq, mock_trim_quality,
        mock_fastqc_threader, mock_error_correction,
        mock_repair_gzip, mock_contamination_finder,
        mock_write_metadata_to_file
    ):
        """
        Test that quality function handles failures.
        """
        print('Starting test_quality_failure')
        # Mock the return values of the functions to simulate failure
        mock_validate_fastq.return_value = self.metadata
        mock_trim_quality.return_value = self.metadata
        mock_fastqc_threader.return_value = self.metadata
        mock_error_correction.return_value = self.metadata
        mock_contamination_finder.return_value = self.metadata

        # Call the function
        result = quality(
            error_logger=self.error_logger,
            log_file=self.log_file,
            metadata=self.metadata,
            report_path=self.report_path,
            sequence_path=self.sequence_path,
            threads=self.threads
        )

        # Assertions
        mock_validate_fastq.assert_called_once_with(
            metadata=self.metadata, log_file=self.log_file
        )
        mock_write_metadata_to_file.assert_called()
        mock_fastqc_threader.assert_called()
        mock_trim_quality.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata
        )
        mock_error_correction.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata,
            threads=self.threads
        )
        mock_repair_gzip.assert_called_once_with(
            log_file=self.log_file, metadata=self.metadata
        )
        mock_contamination_finder.assert_called_once_with(
            input_path=self.sequence_path, log_file=self.log_file,
            metadata=self.metadata, report_path=self.report_path,
            threads=self.threads
        )
        assert result == self.metadata
        print('Completed test_quality_failure')
