# tests/test_pilon.py

"""
Unit tests for pilon misassembly fixing functions in the Teacup COWBAT pipeline
"""

# Standard imports
from unittest.mock import patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.pilon import (
    _execute_pilon_command,
    _initialize_sample_pilon,
    _is_valid_pilon_sample,
    _log_pilon_output,
    _run_pilon,
    pilon
)


class TestPilon:
    """
    Unit tests for Teacup COWBAT pilon functions
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.log_file = 'tests/testdata/log.txt'
        cls.threads = 4
        cls.sample = CustomBox({
            'name': 'sample1',
            'general': {
                'best_assembly_file': 'tests/testdata/sample1_assembly.fa',
                'assembly_file': 'tests/testdata/sample1_assembly.fa',
                'log_out': 'sample1_log_out.txt',
                'log_err': 'sample1_log_err.txt'
            },
            'quast': {
                'sorted_bam': 'tests/testdata/sample1_sorted.bam',
                'output_dir': 'tests/testdata/output/quast'
            }
        })
        cls.metadata = [cls.sample]

    @patch('os.path.isfile')
    @patch('cowbat.pilon._run_pilon')
    @patch('cowbat.pilon._initialize_sample_pilon')
    @patch('cowbat.pilon._is_valid_pilon_sample')
    def test_pilon(
        self, mock_is_valid_pilon_sample, mock_initialize_sample_pilon,
        mock_run_pilon, mock_isfile
    ):
        """
        Test the pilon function.
        """
        mock_is_valid_pilon_sample.return_value = True
        mock_initialize_sample_pilon.return_value = self.sample
        mock_isfile.return_value = False

        updated_metadata = pilon(self.log_file, self.metadata, self.threads)
        mock_is_valid_pilon_sample.assert_called_once_with(self.sample)
        mock_initialize_sample_pilon.assert_called_once_with(
            sample=self.sample, threads=self.threads
        )
        mock_run_pilon.assert_called_once_with(
            log_file=self.log_file, sample=self.sample
        )
        assert updated_metadata == self.metadata

    def test_is_valid_pilon_sample(self):
        """
        Test the _is_valid_pilon_sample function.
        """
        assert _is_valid_pilon_sample(self.sample) is True

        invalid_sample = CustomBox({
            'name': 'sample2',
            'general': {'best_assembly_file': 'NA'}
        })
        assert _is_valid_pilon_sample(invalid_sample) is False

    @patch('os.makedirs')
    def test_initialize_sample_pilon(self, mock_makedirs):
        """
        Test the _initialize_sample_pilon function.
        """
        updated_sample = _initialize_sample_pilon(self.sample, self.threads)
        mock_makedirs.assert_called_once_with(
            'tests/testdata/output/quast/pilon', exist_ok=True
        )
        expected_cmd = (
            'pilon --genome tests/testdata/sample1_assembly.fa '
            '--bam tests/testdata/sample1_sorted.bam --fix bases '
            '--threads 4 --out_dir tests/testdata/output/quast/pilon '
            '--changes --mindepth 0.25'
        )
        assert updated_sample.pilon.cmd == expected_cmd

    @patch('os.path.isfile')
    @patch('cowbat.pilon._execute_pilon_command')
    def test_run_pilon(self, mock_execute_pilon_command, mock_isfile):
        """
        Test the _run_pilon function.
        """
        mock_isfile.return_value = False

        _run_pilon(self.log_file, self.sample)
        mock_isfile.assert_called_once_with(self.sample.general.contigs_file)
        mock_execute_pilon_command.assert_called_once_with(
            log_file=self.log_file, sample=self.sample
        )

    @patch('cowbat.pilon.run_subprocess')
    @patch('cowbat.pilon._log_pilon_output')
    def test_execute_pilon_command(
        self, mock_log_pilon_output, mock_run_subprocess
    ):
        """
        Test the _execute_pilon_command function.
        """
        mock_run_subprocess.return_value = ('output', 'error')

        _execute_pilon_command(self.log_file, self.sample)
        mock_run_subprocess.assert_called_once_with(self.sample.pilon.cmd)
        mock_log_pilon_output.assert_called_once_with(
            err='error',
            log_file=self.log_file,
            out='output',
            sample=self.sample
        )

    @patch('cowbat.pilon.write_to_log_file')
    def test_log_pilon_output(self, mock_write_to_log_file):
        """
        Test the _log_pilon_output function.
        """
        _log_pilon_output(
            err='error',
            log_file=self.log_file,
            out='output', sample=self.sample
        )
        mock_write_to_log_file.assert_called_once_with(
            out=f'{self.sample.pilon.cmd}\noutput', err='error',
            log_file=self.log_file, sample_log=self.sample.general.log_out,
            sample_err=self.sample.general.log_err
        )
