# tests/test_error_correct.py

"""
Unit tests for FASTQ error correction in the Teacup COWBAT pipeline
"""

# Standard imports
from subprocess import CalledProcessError
from unittest.mock import patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
from cowbat.error_correct import error_correction


class TestErrorCorrection:
    """
    Unit tests for Teacup COWBAT error correction functions
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

    @patch('cowbat.error_correct.write_to_log_file')
    @patch('genewrappers.biotools.bbtools.tadpole')
    @patch('os.path.isfile')
    @patch('logging.info')
    @patch('logging.debug')
    def test_error_correction_success(
        self, mock_logging_debug, mock_logging_info,
        mock_isfile, mock_tadpole, mock_write_to_log_file
    ):
        """
        Test that error_correction successfully processes samples.
        """
        print('Starting test_error_correction_success')
        # Mock os.path.isfile to return False (file does not exist)
        mock_isfile.return_value = False

        # Mock bbtools.tadpole to simulate successful execution
        mock_tadpole.return_value = ('stdout', 'stderr', 'cmd')

        # Call the function
        result = error_correction(self.log_file, self.metadata, self.threads)

        # Assertions
        mock_logging_info.assert_any_call('Error correcting reads')
        mock_logging_debug.assert_any_call('Processing sample: %s', 'sample1')
        mock_logging_debug.assert_any_call(
            'Running tadpole for sample: %s', 'sample1'
        )
        mock_logging_debug.assert_any_call('Tadpole command: %s', 'cmd')
        mock_write_to_log_file.assert_called_once_with(
            out='stdout',
            err='stderr',
            log_file=self.log_file,
            sample_log='sample1_log_out.txt',
            sample_err='sample1_log_err.txt'
        )
        assert result[0].general.trimmed_corrected_fastq_files == [
            'sample1_R1_trimmed_corrected.fastq.gz',
            'sample1_R2_trimmed_corrected.fastq.gz'
        ]
        print('Completed test_error_correction_success')

    @patch('os.path.isfile')
    @patch('logging.debug')
    def test_error_correction_file_exists(
        self, mock_logging_debug, mock_isfile
    ):
        """
        Test that error_correction skips samples with existing corrected files.
        """
        print('Starting test_error_correction_file_exists')
        # Mock os.path.isfile to return True (file exists)
        mock_isfile.return_value = True

        # Call the function
        result = error_correction(self.log_file, self.metadata, self.threads)

        # Assertions
        mock_logging_debug.assert_any_call(
            'Trimmed corrected FASTQ file already exists for sample: %s',
            'sample1'
        )
        assert result[0].general.trimmed_corrected_fastq_files == []
        print('Completed test_error_correction_file_exists')

    @patch('genewrappers.biotools.bbtools.tadpole')
    @patch('os.path.isfile')
    @patch('logging.error')
    def test_error_correction_called_process_error(
        self, mock_logging_error, mock_isfile, mock_tadpole
    ):
        """
        Test that error_correction handles CalledProcessError.
        """
        print('Starting test_error_correction_called_process_error')
        # Mock os.path.isfile to return False (file does not exist)
        mock_isfile.return_value = False

        # Mock bbtools.tadpole to raise CalledProcessError
        mock_tadpole.side_effect = CalledProcessError(1, 'cmd')

        # Call the function
        result = error_correction(self.log_file, self.metadata, self.threads)

        # Assertions
        mock_logging_error.assert_any_call(
            'CalledProcessError for sample: %s with error: %s',
            'sample1', 'Command \'cmd\' returned non-zero exit status 1.'
        )
        assert result[0].general.trimmed_corrected_fastq_files == [
            'sample1_R1_trimmed.fastq.gz',
            'sample1_R2_trimmed.fastq.gz'
        ]
        print('Completed test_error_correction_called_process_error')

    @patch('genewrappers.biotools.bbtools.tadpole')
    @patch('os.path.isfile')
    @patch('logging.error')
    def test_error_correction_attribute_error(
        self, mock_logging_error, mock_isfile, mock_tadpole
    ):
        """
        Test that error_correction handles AttributeError.
        """
        print('Starting test_error_correction_attribute_error')
        # Mock os.path.isfile to return False (file does not exist)
        mock_isfile.return_value = False

        # Mock bbtools.tadpole to raise AttributeError
        mock_tadpole.side_effect = AttributeError('attribute error')

        # Call the function
        result = error_correction(self.log_file, self.metadata, self.threads)

        # Assertions
        mock_logging_error.assert_any_call(
            'AttributeError for sample: %s with error: %s',
            'sample1', 'attribute error'
        )
        assert result[0].general.trimmed_corrected_fastq_files == []
        print('Completed test_error_correction_attribute_error')

    @patch('genewrappers.biotools.bbtools.tadpole')
    @patch('os.path.isfile')
    @patch('logging.error')
    def test_error_correction_index_error(
        self, mock_logging_error, mock_isfile, mock_tadpole
    ):
        """
        Test that error_correction handles IndexError.
        """
        print('Starting test_error_correction_index_error')
        # Mock os.path.isfile to return False (file does not exist)
        mock_isfile.return_value = False

        # Mock bbtools.tadpole to raise IndexError
        mock_tadpole.side_effect = IndexError('index error')

        # Call the function
        result = error_correction(self.log_file, self.metadata, self.threads)

        # Assertions
        mock_logging_error.assert_any_call(
            'IndexError for sample: %s with error: %s',
            'sample1', 'index error'
        )
        assert result[0].general.trimmed_corrected_fastq_files == []
        print('Completed test_error_correction_index_error')
