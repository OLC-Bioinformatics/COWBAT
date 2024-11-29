# tests/test_teacup.py

"""
Integration and unit tests for the Teacup COWBAT pipeline
"""

# Standard imports
import os
from unittest.mock import patch, MagicMock

# Third-party imports
import pytest
from olctools.accessoryFunctions.metadata import CustomBox

# Import the Teacup COWBAT pipeline
from cowbat.teacup import TeacupCOWBAT, cli


class TestTeacupCOWBAT:
    """
    Unit tests for the Teacup COWBAT pipeline
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.sequence_path = 'tests/testdata'
        cls.database_path = 'tests/testdata/databases'
        cls.logging_level = 'INFO'
        cls.threads = 4
        cls.teacup = TeacupCOWBAT(
            sequence_path=cls.sequence_path,
            database_path=cls.database_path,
            logging_level=cls.logging_level,
            threads=cls.threads
        )

        # Paths to real FASTQ files for testing
        cls.fastq1 = os.path.join(cls.sequence_path, 'NC_002695_R1.fastq.gz')
        cls.fastq2 = os.path.join(cls.sequence_path, 'NC_002695_R2.fastq.gz')
        cls.json_file = os.path.join(cls.sequence_path, 'NC_002695.json')

        cls.metadata = [
            CustomBox({
                'name': 'sample1',
                'general': {
                    'trimmed_corrected_fastq_files': [
                        cls.fastq1,
                        cls.fastq2
                    ],
                    'fastq_files': [
                        cls.fastq1,
                        cls.fastq2
                    ],
                    'output_directory': os.path.join(
                        cls.sequence_path, 'output'
                    ),
                    'log_out': os.path.join(
                        cls.sequence_path, 'sample1_log_out.txt'
                    ),
                    'log_err': os.path.join(
                        cls.sequence_path, 'sample1_log_err.txt'
                    )
                },
                'commands': CustomBox(),
                'run': CustomBox({
                    'Reads': CustomBox({
                        'forward_read_length': 150,
                        'reverse_read_length': 150
                    })
                }),
                'confindr': CustomBox({
                    'genus': 'Escherichia',
                    'species': 'coli',
                    'contamination_level': 0.01
                }),
                'json_file': cls.json_file
            })
        ]

    def test_initialization(self):
        """
        Test the initialization of the Teacup COWBAT pipeline.
        """
        assert self.teacup.sequence_path == os.path.abspath(
            self.sequence_path
        )
        assert self.teacup.database_path == os.path.abspath(
            self.database_path
        )
        assert self.teacup.report_path == os.path.join(
            os.path.abspath(self.sequence_path), 'reports'
        )
        assert self.teacup.checkpoint_file == os.path.join(
            os.path.abspath(self.sequence_path), '.checkpoint'
        )
        assert self.teacup.log_file == os.path.join(
            os.path.abspath(self.sequence_path), 'log_file'
        )
        assert self.teacup.error_log_file == os.path.join(
            os.path.abspath(self.sequence_path), 'error.log'
        )
        assert self.teacup.threads == self.threads

    @patch('cowbat.teacup.sample_metadata')
    @patch('cowbat.teacup.write_checkpoint')
    @patch('cowbat.teacup.quality', return_value=[])
    @patch('cowbat.teacup.assemble')
    @patch('cowbat.teacup.run_mash_analyses')
    @patch('cowbat.teacup.write_quality_report')
    def test_run_next_step(
        self, mock_write_quality_report, mock_run_mash_analyses,
        mock_assemble, mock_quality, mock_write_checkpoint,
        mock_sample_metadata
    ):
        """
        Test the run_next_step method.
        """
        mock_sample_metadata.return_value = self.metadata
        self.teacup.metadata = mock_sample_metadata()
        self.teacup.last_step = None

        # Debug statement to check metadata before running the step
        print("Metadata before run_next_step:", self.teacup.metadata)

        self.teacup.run_next_step()

        # Debug statement to check metadata after running the step
        print("Metadata after run_next_step:", self.teacup.metadata)

        mock_quality.assert_called_once()
        mock_assemble.assert_called_once()
        mock_run_mash_analyses.assert_called_once()
        mock_write_quality_report.assert_called_once()
        mock_write_checkpoint.assert_called()

    @patch('cowbat.teacup.sample_metadata')
    @patch('cowbat.teacup.write_checkpoint')
    @patch('cowbat.teacup.quality', return_value=[])
    @patch('cowbat.teacup.assemble')
    @patch('cowbat.teacup.run_mash_analyses')
    @patch('cowbat.teacup.write_quality_report')
    def test_main(
        self, mock_write_quality_report, mock_run_mash_analyses,
        mock_assemble, mock_quality, mock_write_checkpoint,
        mock_sample_metadata
    ):
        """
        Test the main method.
        """
        mock_sample_metadata.return_value = self.metadata
        self.teacup.metadata = mock_sample_metadata()

        self.teacup.main()

        mock_sample_metadata.assert_called_once_with(
            error_logger=self.teacup.error_logger,
            metadata=self.teacup.metadata,
            sequence_path=self.teacup.sequence_path
        )
        mock_quality.assert_called_once()
        mock_assemble.assert_called_once()
        mock_run_mash_analyses.assert_called_once()
        mock_write_quality_report.assert_called_once()
        mock_write_checkpoint.assert_called()

    @patch('cowbat.teacup.ArgumentParser.parse_args')
    @patch('cowbat.teacup.TeacupCOWBAT')
    def test_cli(self, mock_teacup_cowbat, mock_parse_args):
        """
        Test the CLI function.
        """
        mock_args = MagicMock()
        mock_args.sequence_path = self.sequence_path
        mock_args.database_path = self.database_path
        mock_args.threads = self.threads
        mock_args.log = self.logging_level
        mock_parse_args.return_value = mock_args

        cli()

        mock_teacup_cowbat.assert_called_once_with(
            sequence_path=self.sequence_path,
            database_path=self.database_path,
            logging_level=self.logging_level,
            threads=self.threads
        )
        mock_teacup_cowbat().main.assert_called_once()


@pytest.mark.integration
def test_teacup_cowbat_integration():
    """
    Integration test for the Teacup COWBAT pipeline.
    """
    sequence_path = 'tests/testdata'
    database_path = 'tests/testdata/databases'
    logging_level = 'INFO'
    threads = 4

    teacup = TeacupCOWBAT(
        sequence_path=sequence_path,
        database_path=database_path,
        logging_level=logging_level,
        threads=threads
    )

    teacup.main()

    # Check if the report file is created
    report_file = os.path.join(sequence_path, 'reports', 'quality_report.tsv')
    assert os.path.exists(report_file)