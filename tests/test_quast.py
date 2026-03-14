# tests/test_quast.py

"""
Unit tests for Quast
"""

# Standard imports
from unittest.mock import patch, mock_open, MagicMock

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.quast import (
    index,
    run_quast,
    prepare_quast_command,
    quast_report_exists,
    run_quast_command,
    parse_quast_report,
    _parse_report,
    clean_quast,
    _remove_large_files,
    analyze
)


class TestQuast:
    """
    Unit tests for assembly evaluation functions.
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.log_file = 'tests/testdata/log.txt'
        cls.threads = 4
        cls.index_queue = MagicMock()
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
            },
            'quast': {
                'sorted_bam': 'tests/testdata/sample1_sorted.bam',
                'output_dir': 'tests/testdata/output/quast',
                'sorted_bai': 'tests/testdata/sample1_sorted.bam.bai',
                'report': 'tests/testdata/output/quast/report.tsv',
                'cmd': (
                    'quast.py --pe1 tests/testdata/sample1_R1.fastq.gz '
                    '--pe2 tests/testdata/sample1_R2.fastq.gz '
                    '--ref-bam tests/testdata/sample1_sorted.bam -t 4 '
                    '--k-mer-stats --circos --rna-finding '
                    '--conserved-genes-finding -o tests/testdata/output/quast '
                    '--debug tests/testdata/sample1_assembly.fa --threads 4'
                )
            }
        })
        cls.metadata = [cls.sample]

    @patch('os.path.isfile')
    def test_index(self, mock_isfile):
        """
        Test the index function.
        """
        mock_isfile.return_value = False
        bam_index = MagicMock()
        bam_index.return_value = ("stdout", "stderr")
        self.index_queue.get.side_effect = [
            (self.sample, bam_index), None
        ]
        index(self.index_queue)
        self.index_queue.get.assert_called()
        self.index_queue.task_done.assert_called()

    @patch('os.path.isfile')
    def test_quast_report_exists(self, mock_isfile):
        """
        Test the quast_report_exists function.
        """
        mock_isfile.return_value = True
        assert quast_report_exists(self.sample) is True

    @patch('cowbat.quast.run_subprocess')
    @patch('cowbat.quast.write_to_log_file')
    def test_run_quast_command(
        self, mock_write_to_log_file, mock_run_subprocess
    ):
        """
        Test the run_quast_command function.
        """
        mock_run_subprocess.return_value = ('output', 'error')
        run_quast_command(self.sample, self.log_file)
        mock_run_subprocess.assert_called_once_with(self.sample.quast.cmd)
        mock_write_to_log_file.assert_called_once()

    @patch('cowbat.quast.prepare_quast_command')
    @patch('cowbat.quast.quast_report_exists')
    @patch('cowbat.quast.run_quast_command')
    def test_run_quast(
        self, mock_run_quast_command, mock_quast_report_exists,
        mock_prepare_quast_command
    ):
        """
        Test the run_quast function.
        """
        mock_quast_report_exists.return_value = False

        # Ensure it returns the expected CustomBox object
        mock_prepare_quast_command.return_value = self.sample

        run_quast(self.log_file, self.metadata, self.threads)
        mock_prepare_quast_command.assert_called_once_with(
            sample=self.sample, threads=self.threads
        )
        mock_quast_report_exists.assert_called_once_with(self.sample)
        mock_run_quast_command.assert_called_once_with(
            sample=self.sample, log_file=self.log_file
        )

    def test_prepare_quast_command(self):
        """
        Test the prepare_quast_command function.
        """
        updated_sample = prepare_quast_command(self.sample, self.threads)
        expected_command = (
            'quast.py --pe1 tests/testdata/sample1_R1.fastq.gz '
            '--pe2 tests/testdata/sample1_R2.fastq.gz '
            '--ref-bam tests/testdata/sample1_sorted.bam -t 4 '
            '--k-mer-stats --circos --rna-finding '
            '--conserved-genes-finding -o tests/testdata/output/quast '
            '--debug tests/testdata/sample1_assembly.fa --threads 4'
        )
        assert updated_sample.quast.cmd == expected_command

    @patch('os.path.isfile')
    def test_parse_quast_report(self, mock_isfile):
        """
        Test the parse_quast_report function.
        """
        mock_isfile.return_value = True
        updated_metadata = parse_quast_report(self.metadata)
        assert updated_metadata == self.metadata

    @patch('cowbat.quast.analyze')
    @patch('builtins.open', new_callable=mock_open, read_data='key\tvalue\n')
    def test__parse_report(self, _, mock_analyze):
        """
        Test the _parse_report function.
        """
        mock_analyze.return_value = ('key', 'value')
        updated_sample = _parse_report(self.sample)
        assert getattr(updated_sample.quast, 'key') == 'value'

    @patch('os.path.isdir')
    @patch('cowbat.quast._clean_sample_quast_files')
    def test_clean_quast(
        self, mock_clean_sample_quast_files, mock_isdir
    ):
        """
        Test the clean_quast function.
        """
        mock_isdir.return_value = True
        clean_quast(self.metadata)
        mock_clean_sample_quast_files.assert_called_once_with(self.sample)

    @patch('os.path.getsize')
    @patch('os.remove')
    def test__remove_large_files(self, mock_remove, mock_getsize):
        """
        Test the _remove_large_files function.
        """
        mock_getsize.return_value = 200000
        _remove_large_files(
            'tests/testdata/large_file.txt', 'large_file.txt', self.sample
        )
        mock_remove.assert_called_once_with('tests/testdata/large_file.txt')

    def test_analyze(self):
        """
        Test the analyze function.
        """
        line = "key\tvalue"
        key, value = analyze(line)
        assert key == "key"
        assert value == "value"
