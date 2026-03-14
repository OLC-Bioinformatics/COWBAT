# tests/test_quality_trim.py

"""
Unit tests for FASTQ quality trimming in the Teacup COWBAT pipeline
"""

# Standard imports
from subprocess import CalledProcessError
from unittest.mock import patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
from cowbat.quality_trim import trim_quality


class TestQualityTrim:
    """
    Unit tests for Teacup COWBAT quality trimming functions
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
                    'fastq_files': [
                        'sample1_R1.fastq.gz', 'sample1_R2.fastq.gz'
                    ],
                    'trimmed_fastq_files': [],
                    'log_out': 'sample1_log_out.txt',
                    'log_err': 'sample1_log_err.txt',
                    'output_directory': 'tests/testdata/output'
                },
                'run': {
                    'Reads': {
                        'forward_read_length': '150',
                        'reverse_read_length': '150'
                    }
                },
                'commands': {}
            })
        ]

    @patch('cowbat.quality_trim.bbtools.bbduk_trim')
    @patch('cowbat.quality_trim.write_to_log_file')
    @patch('cowbat.quality_trim.glob')
    @patch('os.path.isfile')
    def test_trim_quality_success(
        self, mock_isfile, mock_glob, mock_write_to_log_file, mock_bbduk_trim
    ):
        """
        Test that trim_quality successfully processes samples.
        """
        print('Starting test_trim_quality_success')
        mock_isfile.return_value = False
        mock_glob.return_value = [
            'tests/testdata/output/sample1_R1_trimmed.fastq.gz',
            'tests/testdata/output/sample1_R2_trimmed.fastq.gz'
        ]
        mock_bbduk_trim.return_value = ('stdout', 'stderr', 'bbduk_call')

        result = trim_quality(self.log_file, self.metadata)

        mock_bbduk_trim.assert_called()
        mock_write_to_log_file.assert_called()
        assert result[0].general.trimmed_fastq_files == mock_glob.return_value
        print('Completed test_trim_quality_success')

    @patch('os.path.isfile')
    def test_trim_quality_file_exists(self, mock_isfile):
        """
        Test that trim_quality skips samples with existing trimmed files.
        """
        print('Starting test_trim_quality_file_exists')
        mock_isfile.return_value = True

        result = trim_quality(self.log_file, self.metadata)

        assert result[0].general.trimmed_fastq_files == []
        print('Completed test_trim_quality_file_exists')

    @patch('cowbat.quality_trim.bbtools.bbduk_trim')
    @patch('os.path.isfile')
    def test_trim_quality_called_process_error(
        self, mock_isfile, mock_bbduk_trim
    ):
        """
        Test that trim_quality handles CalledProcessError.
        """
        print('Starting test_trim_quality_called_process_error')
        mock_isfile.return_value = False
        mock_bbduk_trim.side_effect = CalledProcessError(1, 'cmd')

        result = trim_quality(self.log_file, self.metadata)

        assert result[0].general.trimmed_fastq_files == []
        print('Completed test_trim_quality_called_process_error')

    @patch('cowbat.quality_trim.bbtools.bbduk_trim')
    @patch('os.path.isfile')
    def test_trim_quality_attribute_error(self, mock_isfile, mock_bbduk_trim):
        """
        Test that trim_quality handles AttributeError.
        """
        print('Starting test_trim_quality_attribute_error')
        mock_isfile.return_value = False
        mock_bbduk_trim.side_effect = AttributeError('attribute error')

        result = trim_quality(self.log_file, self.metadata)

        assert result[0].general.trimmed_fastq_files == []
        print('Completed test_trim_quality_attribute_error')

    @patch('cowbat.quality_trim.bbtools.bbduk_trim')
    @patch('os.path.isfile')
    def test_trim_quality_index_error(self, mock_isfile, mock_bbduk_trim):
        """
        Test that trim_quality handles IndexError.
        """
        print('Starting test_trim_quality_index_error')
        mock_isfile.return_value = False
        mock_bbduk_trim.side_effect = IndexError('index error')

        result = trim_quality(self.log_file, self.metadata)

        assert result[0].general.trimmed_fastq_files == []
        print('Completed test_trim_quality_index_error')
