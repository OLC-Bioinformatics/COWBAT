# tests/test_prodigal.py

"""
Unit tests for Prodigal gene prediction functions in the Teacup COWBAT pipeline
"""

# Standard imports
from unittest.mock import patch, MagicMock

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Import the Prodigal class from cowbat.prodigal
from cowbat.prodigal import Prodigal


class TestProdigal:
    """
    Unit tests for Teacup COWBAT Prodigal functions
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
                    'best_assembly_file': 'tests/testdata/sample1_assembly.fa',
                    'output_directory': 'tests/testdata/output',
                    'log_out': 'sample1_log_out.txt',
                    'log_err': 'sample1_log_err.txt'
                },
                'commands': {},
                'prodigal': CustomBox()
            })
        ]
        cls.prodigal = Prodigal(cls.log_file, cls.metadata)

    def test_populate_prodigal_attributes(self):
        """
        Test the populate_prodigal_attributes function.
        """
        sample = self.metadata[0]
        self.prodigal.populate_prodigal_attributes(sample)
        assert sample.prodigal.report_dir == 'tests/testdata/output/prodigal'
        assert sample.prodigal.results_file == (
            'tests/testdata/output/prodigal/sample1_prodigalresults.sco'
        )
        assert sample.commands.prodigal == (
            'prodigal -i tests/testdata/sample1_assembly.fa '
            '-o tests/testdata/output/prodigal/sample1_prodigalresults.sco '
            '-f sco -d tests/testdata/output/prodigal/sample1_genes.fa'
        )

    @patch('os.makedirs')
    def test_create_report_directory(self, mock_makedirs):
        """
        Test the create_report_directory function.
        """
        sample = self.metadata[0]
        self.prodigal.create_report_directory(sample)
        mock_makedirs.assert_called_once_with(
            'tests/testdata/output/prodigal', exist_ok=True
        )

    @patch('os.path.isfile')
    @patch('os.stat')
    def test_is_report_needed(self, mock_stat, mock_isfile):
        """
        Test the is_report_needed function.
        """
        sample = self.metadata[0]
        sample.prodigal.results_file = (
            'tests/testdata/sample1_prodigalresults.sco'
        )

        # Case 1: File does not exist
        mock_isfile.return_value = False
        assert self.prodigal.is_report_needed(sample) is True

        # Case 2: File exists but is empty
        mock_isfile.return_value = True
        mock_stat.return_value.st_size = 0
        assert self.prodigal.is_report_needed(sample) is True

        # Case 3: File exists and is not empty
        mock_stat.return_value.st_size = 100
        assert self.prodigal.is_report_needed(sample) is False

    @patch('cowbat.prodigal.run_subprocess')
    @patch('cowbat.prodigal.write_to_log_file')
    def test_run_prodigal_command(
        self, mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test the run_prodigal_command function.
        """
        sample = self.metadata[0]
        sample.commands.prodigal = 'prodigal command'
        mock_run_subprocess.return_value = ('output', 'error')

        self.prodigal.run_prodigal_command(sample)
        mock_run_subprocess.assert_called_once_with('prodigal command')
        mock_write_to_log_file.assert_any_call(
            out='prodigal command',
            err='prodigal command',
            log_file=self.log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )
        mock_write_to_log_file.assert_any_call(
            out='output',
            err='error',
            log_file=self.log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )

    def test_initialize_prodigal_attributes(self):
        """
        Test the initialize_prodigal_attributes function.
        """
        sample = self.metadata[0]
        self.prodigal.initialize_prodigal_attributes(sample)
        assert sample.prodigal.predicted_genes_total == 0
        assert sample.prodigal.predicted_genes_over_3000bp == 0
        assert sample.prodigal.predicted_genes_over_1000bp == 0
        assert sample.prodigal.predicted_genes_over_500bp == 0
        assert sample.prodigal.predicted_genes_under_500bp == 0

    @patch('builtins.open', new_callable=MagicMock)
    def test_parse_prodigal_results(self, mock_open):
        """
        Test the parse_prodigal_results function.
        """
        sample = self.metadata[0]
        sample.prodigal.results = (
            'tests/testdata/sample1_prodigalresults.sco'
        )
        mock_open.return_value.__enter__.return_value = [
            '>gene_1_100_200', '>gene_2_300_350', '>gene_3_400_450'
        ]

        self.prodigal.parse_prodigal_results(sample)
        assert sample.prodigal.predicted_genes_total == 3
        assert sample.prodigal.predicted_genes_over_3000bp == 0
        assert sample.prodigal.predicted_genes_over_1000bp == 0
        assert sample.prodigal.predicted_genes_over_500bp == 0
        assert sample.prodigal.predicted_genes_under_500bp == 3
