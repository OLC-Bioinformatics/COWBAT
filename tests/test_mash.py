# tests/test_mash.py

"""
Unit tests for Mash analysis functions in the Teacup COWBAT pipeline
"""

# Standard imports
import logging
import os
from unittest.mock import mock_open, patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Import functions directly from genemethods.assemblypipeline.mash
from genemethods.assemblypipeline.mash import (
    _check_file_extension,
    _create_fastq_list_file,
    _create_metadata_attributes,
    _load_mash_data,
    _parse_mash_data,
    _populate_ref_dict,
    _run_mash_sketch,
    mashing,
    parsing,
    run_mash_analyses,
    sketching,
    write_mash_report
)


class TestMash:
    """
    Unit tests for Teacup COWBAT Mash functions
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.analysis_type = 'mash_analysis'
        cls.log_file = 'tests/testdata/log.txt'
        cls.reference_file_path = 'tests/testdata/reference'
        cls.report_path = 'tests/testdata/report'
        cls.threads = 4
        cls.metadata = [
            CustomBox({
                'name': 'sample1',
                'general': {
                    'output_directory': 'tests/testdata/output',
                    'trimmed_corrected_fastq_files': [
                        'tests/testdata/sample1.fastq'
                    ],
                    'log_out': 'sample1_log_out.txt',
                    'log_err': 'sample1_log_err.txt'
                },
                'commands': {},
                'mash_analysis': CustomBox()
            })
        ]

    @patch('genemethods.assemblypipeline.mash.sketching')
    @patch('genemethods.assemblypipeline.mash.mashing')
    @patch('genemethods.assemblypipeline.mash.parsing')
    @patch('genemethods.assemblypipeline.mash.write_mash_report')
    @patch('genemethods.assemblypipeline.mash.write_metadata_to_file')
    def test_run_mash_analyses(
        self, mock_write_metadata_to_file, mock_write_mash_report,
        mock_parsing, mock_mashing, mock_sketching
    ):
        """
        Test the run_mash_analyses function.
        """
        mock_sketching.return_value = self.metadata
        mock_parsing.return_value = self.metadata

        updated_metadata = run_mash_analyses(
            analysis_type=self.analysis_type,
            error_logger=logging,
            log_file=self.log_file,
            metadata=self.metadata,
            reference_file_path=self.reference_file_path,
            report_path=self.report_path,
            threads=self.threads
        )

        mock_sketching.assert_called_once_with(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            metadata=self.metadata,
            reference_file_path=os.path.join(
                self.reference_file_path, self.analysis_type
            ),
            threads=self.threads
        )
        mock_mashing.assert_called_once_with(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            metadata=self.metadata
        )
        mock_parsing.assert_called_once_with(
            analysis_type=self.analysis_type,
            metadata=self.metadata,
            reference_file_path=os.path.join(
                self.reference_file_path, self.analysis_type
            )
        )
        mock_write_mash_report.assert_called_once_with(
            analysis_type=self.analysis_type,
            metadata=self.metadata,
            report_path=self.report_path
        )
        mock_write_metadata_to_file.assert_called_once_with(
            error_logger=logging,
            metadata=self.metadata
        )
        assert updated_metadata == self.metadata

    @patch('genemethods.assemblypipeline.mash._create_metadata_attributes')
    @patch('genemethods.assemblypipeline.mash._create_fastq_list_file')
    @patch('genemethods.assemblypipeline.mash._run_mash_sketch')
    def test_sketching(
        self, mock_run_mash_sketch, mock_create_fastq_list_file,
        mock_create_metadata_attributes
    ):
        """
        Test the sketching function.
        """
        mock_create_metadata_attributes.return_value = self.metadata[0]

        updated_metadata = sketching(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            metadata=self.metadata,
            reference_file_path=self.reference_file_path,
            threads=self.threads
        )

        mock_create_metadata_attributes.assert_called_once_with(
            analysis_type=self.analysis_type,
            reference_file_path=self.reference_file_path,
            sample=self.metadata[0],
            threads=self.threads
        )
        mock_create_fastq_list_file.assert_called_once_with(
            analysis_type=self.analysis_type,
            sample=self.metadata[0]
        )
        mock_run_mash_sketch.assert_called_once_with(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            sample=self.metadata[0]
        )
        assert updated_metadata == self.metadata

    def test_create_metadata_attributes(self):
        """
        Test the _create_metadata_attributes function.
        """
        sample = self.metadata[0]
        updated_sample = _create_metadata_attributes(
            analysis_type=self.analysis_type,
            reference_file_path=self.reference_file_path,
            sample=sample,
            threads=self.threads
        )
        assert updated_sample[self.analysis_type].report_dir == (
            'tests/testdata/output/mash_analysis'
        )
        assert updated_sample[self.analysis_type].sketch_file == (
            'tests/testdata/output/mash_analysis/sample1.msh'
        )
        assert updated_sample.commands.sketch.startswith('mash sketch')
        assert updated_sample.commands.mash.startswith('mash dist')

    @patch('builtins.open', new_callable=mock_open)
    def test_create_fastq_list_file(self, mock_open_instance):
        """
        Test the _create_fastq_list_file function.
        """
        sample = self.metadata[0]
        _create_fastq_list_file(
            analysis_type=self.analysis_type,
            sample=sample
        )
        mock_open_instance.assert_called_once_with(
            sample[self.analysis_type].file_list, 'w', encoding='utf-8'
        )
        mock_open_instance().write.assert_called_once_with(
            'tests/testdata/sample1.fastq'
        )

    def test_check_file_extension(self):
        """
        Test the _check_file_extension function.
        """
        sample = self.metadata[0]
        assert _check_file_extension(sample) is False
        sample.general.trimmed_corrected_fastq_files = [
            'tests/testdata/sample1.fasta'
        ]
        assert _check_file_extension(sample) is True

    @patch('os.path.isfile')
    @patch('genemethods.assemblypipeline.mash.run_subprocess')
    @patch('genemethods.assemblypipeline.mash.write_to_log_file')
    def test_run_mash_sketch(
        self, mock_write_to_log_file, mock_run_subprocess, mock_isfile
    ):
        """
        Test the _run_mash_sketch function.
        """
        sample = self.metadata[0]
        sample.commands.sketch = 'mash sketch command'
        mock_isfile.return_value = False
        mock_run_subprocess.return_value = ('output', 'error')

        _run_mash_sketch(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            sample=sample
        )
        mock_run_subprocess.assert_called_once_with('mash sketch command')
        mock_write_to_log_file.assert_any_call(
            out='mash sketch command',
            err='mash sketch command',
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

    @patch('os.path.isfile')
    @patch('genemethods.assemblypipeline.mash.run_subprocess')
    @patch('genemethods.assemblypipeline.mash.write_to_log_file')
    def test_mashing(
        self, mock_write_to_log_file, mock_run_subprocess, mock_isfile
    ):
        """
        Test the mashing function.
        """
        sample = self.metadata[0]
        sample.commands.mash = 'mash command'
        mock_isfile.return_value = False
        mock_run_subprocess.return_value = ('output', 'error')

        mashing(
            analysis_type=self.analysis_type,
            log_file=self.log_file,
            metadata=self.metadata
        )
        mock_run_subprocess.assert_called_once_with('mash command')
        mock_write_to_log_file.assert_any_call(
            out='mash command',
            err='mash command',
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

    @patch('builtins.open', new_callable=mock_open, read_data='data')
    def test_populate_ref_dict(self, mock_open_instance):
        """
        Test the _populate_ref_dict function.
        """
        ref_dict = _populate_ref_dict(self.reference_file_path)
        mock_open_instance.assert_called_once_with(
            os.path.join(
                self.reference_file_path, 'assembly_summary_refseq.txt'
            ),
            encoding='utf-8'
        )
        assert isinstance(ref_dict, dict)

    @patch('builtins.open', new_callable=mock_open, read_data='data')
    def test_load_mash_data(self, mock_open_instance):
        """
        Test the _load_mash_data function.
        """
        sample = self.metadata[0]
        sample[self.analysis_type].mash_results = 'tests/testdata/sample1.tab'
        mash_data = _load_mash_data(
            analysis_type=self.analysis_type,
            sample=sample
        )
        mock_open_instance.assert_called_once_with(
            'tests/testdata/sample1.tab', 'r', encoding='utf-8'
        )
        assert isinstance(mash_data, list)

    @patch('genemethods.assemblypipeline.mash._populate_ref_dict')
    @patch('genemethods.assemblypipeline.mash._load_mash_data')
    @patch('genemethods.assemblypipeline.mash._parse_mash_data')
    def test_parsing(
        self, mock_parse_mash_data, mock_load_mash_data,
        mock_populate_ref_dict
    ):
        """
        Test the parsing function.
        """
        mock_populate_ref_dict.return_value = {
            "GCF_000008865": "Helicobacter pullorum"
        }
        mock_load_mash_data.return_value = [
            "GCF_000008865.1\t0.001\t0.05\t100"
        ]
        mock_parse_mash_data.return_value = self.metadata[0]

        updated_metadata = parsing(
            analysis_type=self.analysis_type,
            metadata=self.metadata,
            reference_file_path=self.reference_file_path
        )

        mock_populate_ref_dict.assert_called_once_with(
            reference_file_path=self.reference_file_path
        )
        mock_load_mash_data.assert_called_once_with(
            analysis_type=self.analysis_type,
            sample=self.metadata[0]
        )
        mock_parse_mash_data.assert_called_once_with(
            analysis_type=self.analysis_type,
            mash_data=["GCF_000008865.1\t0.001\t0.05\t100"],
            ref_dict={"GCF_000008865": "Helicobacter pullorum"},
            sample=self.metadata[0]
        )
        assert updated_metadata == self.metadata

    def test_parse_mash_data(self):
        """
        Test the _parse_mash_data function.
        """
        sample = self.metadata[0]
        mash_data = ["GCF_000008865.1\t000\t0.001\t0.05\t100"]
        ref_dict = {"GCF_000008865": "Helicobacter pullorum"}
        updated_sample = _parse_mash_data(
            analysis_type=self.analysis_type,
            mash_data=mash_data,
            ref_dict=ref_dict,
            sample=sample
        )
        assert updated_sample[self.analysis_type].closest_refseq == (
            "Helicobacter pullorum"
        )
        assert updated_sample[self.analysis_type].closest_refseq_genus == (
            "Helicobacter"
        )
        assert updated_sample[self.analysis_type].closest_refseq_species == (
            "pullorum"
        )
        assert updated_sample[self.analysis_type].mash_distance == "0.001"
        assert updated_sample[self.analysis_type].p_value == "0.05"
        assert updated_sample[self.analysis_type].num_matches == "100"

    @patch('builtins.open', new_callable=mock_open)
    def test_write_mash_report(self, mock_open_instance):
        """
        Test the write_mash_report function.
        """
        # Update the metadata to reflect the expected values
        self.metadata[0][self.analysis_type].closest_refseq = \
            "Helicobacter pullorum"
        self.metadata[0][self.analysis_type].closest_refseq_genus = \
            "Helicobacter"
        self.metadata[0][self.analysis_type].closest_refseq_species = \
            "pullorum"
        self.metadata[0][self.analysis_type].mash_distance = "0.001"
        self.metadata[0][self.analysis_type].p_value = "0.05"
        self.metadata[0][self.analysis_type].num_matches = "100"

        write_mash_report(
            analysis_type=self.analysis_type,
            metadata=self.metadata,
            report_path=self.report_path
        )
        mock_open_instance.assert_called_once_with(
            os.path.join(self.report_path, 'mash.tsv'), 'w', encoding='utf-8'
        )
        mock_open_instance().write.assert_any_call(
            'Strain\treference_genus\treference_file\t'
            'reference_genome_mash_distance\tp_value\tnum_matching_hashes\n'
        )
        mock_open_instance().write.assert_any_call(
            'sample1\tHelicobacter\tHelicobacter pullorum\t0.001\t0.05\t100\n'
        )
