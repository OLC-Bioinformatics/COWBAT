# tests/test_bowtie2.py

"""
Unit tests for Bowtie2 and Samtools mapping functions in the Teacup COWBAT
pipeline.
"""

# Standard imports
from unittest.mock import patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.bowtie2 import (
    prepare_sample_for_build,
    index_files_exist,
    run_build_command,
    bowtie_build,
    prepare_output_directory,
    construct_mapping_command,
    sorted_bam_exists,
    run_mapping_command,
    bowtie_run
)


class TestBowtie2:
    """
    Unit tests for Bowtie2 and Samtools mapping functions.
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
                'assembly_file': 'tests/testdata/sample1_assembly.fa',
                'trimmed_corrected_fastq_files': [
                    'tests/testdata/sample1_R1.fastq.gz',
                    'tests/testdata/sample1_R2.fastq.gz'
                ],
                'output_directory': 'tests/testdata/output',
                'log_out': 'sample1_log_out.txt',
                'log_err': 'sample1_log_err.txt',
                'best_assembly_file': 'tests/testdata/sample1_assembly.fa'
            }
        })
        cls.metadata = [cls.sample]

    def test_prepare_sample_for_build(self):
        """
        Test the prepare_sample_for_build function.
        """
        prepare_sample_for_build(self.sample)
        assert self.sample.quast.base_name == 'tests/testdata/sample1_assembly'
        assert self.sample.quast.build_command == (
            'bowtie2-build tests/testdata/sample1_assembly.fa '
            'tests/testdata/sample1_assembly'
        )

    @patch('os.path.isfile')
    def test_index_files_exist(self, mock_isfile):
        """
        Test the index_files_exist function.
        """
        mock_isfile.return_value = True
        prepare_sample_for_build(self.sample)
        assert index_files_exist(self.sample) is True

    @patch('cowbat.bowtie2.run_subprocess')
    @patch('cowbat.bowtie2.write_to_log_file')
    def test_run_build_command(
        self, mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test the run_build_command function.
        """
        mock_run_subprocess.return_value = ('output', 'error')
        prepare_sample_for_build(self.sample)
        run_build_command(self.sample, self.log_file)
        mock_run_subprocess.assert_called_once_with(
            self.sample.quast.build_command
        )
        mock_write_to_log_file.assert_called_once()

    @patch('cowbat.bowtie2.run_build_command')
    @patch('cowbat.bowtie2.index_files_exist')
    @patch('cowbat.bowtie2.prepare_sample_for_build')
    def test_bowtie_build(
        self, mock_prepare_sample_for_build, mock_index_files_exist,
        mock_run_build_command
    ):
        """
        Test the bowtie_build function.
        """
        mock_index_files_exist.return_value = False
        bowtie_build(self.log_file, self.metadata)
        mock_prepare_sample_for_build.assert_called_once_with(self.sample)
        mock_index_files_exist.assert_called_once_with(self.sample)
        mock_run_build_command.assert_called_once_with(
            self.sample, self.log_file
        )

    def test_prepare_output_directory(self):
        """
        Test the prepare_output_directory function.
        """
        prepare_output_directory(self.sample)
        assert self.sample.quast.output_dir == 'tests/testdata/output/quast'
        assert self.sample.quast.sorted_bam == (
            'tests/testdata/output/quast/sample1_sorted.bam'
        )

    def test_construct_mapping_command(self):
        """
        Test the construct_mapping_command function.
        """
        prepare_sample_for_build(self.sample)
        prepare_output_directory(self.sample)
        construct_mapping_command(self.sample, self.threads)
        expected_command = (
            'bowtie2 -x tests/testdata/sample1_assembly '
            '-1 tests/testdata/sample1_R1.fastq.gz '
            '-2 tests/testdata/sample1_R2.fastq.gz '
            '-p 4 -X 1000 | '
            'samtools view -@ 4 -h -F 4 -bT '
            'tests/testdata/sample1_assembly.fa - | '
            'samtools sort - -@ 4 -o '
            'tests/testdata/output/quast/sample1_sorted.bam'
        )
        assert self.sample.quast.map_command == expected_command

    @patch('os.path.isfile')
    def test_sorted_bam_exists(self, mock_isfile):
        """
        Test the sorted_bam_exists function.
        """
        mock_isfile.return_value = True
        prepare_output_directory(self.sample)
        assert sorted_bam_exists(self.sample) is True

    @patch('cowbat.bowtie2.run_subprocess')
    @patch('cowbat.bowtie2.write_to_log_file')
    def test_run_mapping_command(
        self, mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test the run_mapping_command function.
        """
        mock_run_subprocess.return_value = ('output', 'error')
        prepare_sample_for_build(self.sample)
        prepare_output_directory(self.sample)
        construct_mapping_command(self.sample, self.threads)
        run_mapping_command(self.sample, self.log_file)
        mock_run_subprocess.assert_called_once_with(
            self.sample.quast.map_command
        )
        mock_write_to_log_file.assert_called_once()

    @patch('cowbat.bowtie2.run_mapping_command')
    @patch('cowbat.bowtie2.sorted_bam_exists')
    @patch('cowbat.bowtie2.construct_mapping_command')
    @patch('cowbat.bowtie2.prepare_output_directory')
    def test_bowtie_run(
        self, mock_prepare_output_directory, mock_construct_mapping_command,
        mock_sorted_bam_exists, mock_run_mapping_command
    ):
        """
        Test the bowtie_run function.
        """
        mock_sorted_bam_exists.return_value = False
        bowtie_run(self.log_file, self.metadata, self.threads)
        mock_prepare_output_directory.assert_called_once_with(self.sample)
        mock_construct_mapping_command.assert_called_once_with(
            self.sample, self.threads
        )
        mock_sorted_bam_exists.assert_called_once_with(self.sample)
        mock_run_mapping_command.assert_called_once_with(
            self.sample, self.log_file
        )
