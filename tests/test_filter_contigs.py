# tests/test_filter_contigs.py

"""
Unit tests for contig filtering functions in the Teacup COWBAT pipeline
"""

# Standard imports
import os
from unittest.mock import patch, mock_open, MagicMock

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Import functions directly from cowbat.filter
from cowbat.filter import (
    _copy_filtered_file,
    _filter_contigs,
    _initialize_sample_filter,
    _is_valid_filter_sample,
    _run_filter
)


class TestFilterContigs:
    """
    Unit tests for Teacup COWBAT contig filtering functions
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.sequence_path = 'tests/testdata/sequence'
        cls.metadata = CustomBox({
            'name': 'sample1',
            'general': {
                'best_assembly_file': 'tests/testdata/sample1_assembly.fa',
                'assembly_file': 'tests/testdata/sample1_assembly.fa',
                'filtered_file': 'tests/testdata/sample1_filtered.fa'
            },
            'qualimap': {
                'coverage': {'contig1': '30.0'},
                'std_dev': {'contig1': '5.0'}
            }
        })

    def test_is_valid_filter_sample(self):
        """
        Test the _is_valid_filter_sample function.
        """
        assert _is_valid_filter_sample(self.metadata) is True

        invalid_sample = CustomBox({
            'name': 'sample2',
            'general': {'best_assembly_file': 'NA'}
        })
        assert _is_valid_filter_sample(invalid_sample) is False

    def test_initialize_sample_filter(self):
        """
        Test the _initialize_sample_filter function.
        """
        updated_sample = _initialize_sample_filter(self.metadata)
        assert updated_sample.general.contigs_file == (
            self.metadata.general.assembly_file
        )

    @patch('os.path.isfile')
    @patch('cowbat.filter._copy_filtered_file')
    @patch('cowbat.filter._filter_contigs')
    def test_run_filter(
        self, mock_filter_contigs, mock_copy_filtered_file, mock_isfile
    ):
        """
        Test the _run_filter function.
        """
        mock_isfile.return_value = True
        mock_copy_filtered_file.return_value = self.metadata

        updated_sample = _run_filter(self.metadata, self.sequence_path)
        mock_isfile.assert_called_once_with(self.metadata.general.contigs_file)
        mock_filter_contigs.assert_called_once_with(self.metadata)
        mock_copy_filtered_file.assert_called_once_with(
            sample=self.metadata, sequence_path=self.sequence_path
        )
        assert updated_sample == self.metadata

    @patch(
        'builtins.open', new_callable=mock_open, read_data='>contig1\nATGC\n'
    )
    @patch('Bio.SeqIO.parse')
    def test_filter_contigs(self, mock_seqio_parse, mock_open_instance):
        """
        Test the _filter_contigs function.
        """
        mock_seqio_parse.return_value = [
            MagicMock(id='contig1_pilon', seq='ATGC')
        ]

        _filter_contigs(self.metadata)
        mock_open_instance.assert_called_once_with(
            self.metadata.general.contigs_file, 'r', encoding='utf-8'
        )
        mock_seqio_parse.assert_called_once_with(mock_open_instance(), 'fasta')

    @patch('os.path.isfile')
    @patch('shutil.copyfile')
    def test_copy_filtered_file(self, mock_copyfile, mock_isfile):
        """
        Test the _copy_filtered_file function.
        """
        mock_isfile.side_effect = [True, False]

        updated_sample = _copy_filtered_file(self.metadata, self.sequence_path)
        mock_isfile.assert_any_call(self.metadata.general.filtered_file)
        mock_isfile.assert_any_call(self.metadata.general.best_assembly_file)
        mock_copyfile.assert_called_once_with(
            self.metadata.general.filtered_file,
            self.metadata.general.best_assembly_file
        )
        assert updated_sample.general.best_assembly_file == os.path.join(
            self.sequence_path, 'BestAssemblies', f'{self.metadata.name}.fasta'
        )
