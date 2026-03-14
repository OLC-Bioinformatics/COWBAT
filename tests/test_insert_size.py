# tests/test_insert_size.py

"""
Unit tests for insert size estimation functions in the Teacup COWBAT pipeline
"""

# Standard imports
from unittest.mock import patch, mock_open

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox
import cowbat.insert_size as insert_size


class TestInsertSize:
    """
    Unit tests for Teacup COWBAT insert size functions
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.sample = CustomBox({
            'name': 'sample1',
            'quast': {
                'output_dir': 'tests/testdata/output/quast',
                'reads_stats_file': (
                    'tests/testdata/output/quast/reads_stats/reads_stats.err'
                ),
                'total_reads': 0
            }
        })
        cls.metadata = [cls.sample]

    @patch('os.path.isfile')
    @patch('cowbat.insert_size.initialize_sample_attributes')
    @patch('cowbat.insert_size.parse_reads_stats_file')
    def test_extract_insert_size(
        self, mock_parse_reads_stats_file, mock_initialize_sample_attributes,
        mock_isfile
    ):
        """
        Test the extract_insert_size function.
        """
        mock_isfile.return_value = True
        mock_initialize_sample_attributes.return_value = self.sample
        mock_parse_reads_stats_file.return_value = self.sample

        updated_metadata = insert_size.extract_insert_size(self.metadata)
        mock_isfile.assert_called_once_with(self.sample.quast.reads_stats_file)
        mock_initialize_sample_attributes.assert_called_once_with(self.sample)
        mock_parse_reads_stats_file.assert_called_once_with(self.sample)
        assert updated_metadata == self.metadata

    @patch('os.path.isfile')
    def test_extract_insert_size_file_not_found(self, mock_isfile):
        """
        Test the extract_insert_size function when the reads stats file is not
        found.
        """
        mock_isfile.return_value = False

        updated_metadata = insert_size.extract_insert_size(self.metadata)
        mock_isfile.assert_called_once_with(self.sample.quast.reads_stats_file)
        assert self.sample.quast.insert_mean == 'ND'
        assert self.sample.quast.insert_std == 'ND'
        assert updated_metadata == self.metadata

    def test_initialize_sample_attributes(self):
        """
        Test the initialize_sample_attributes function.
        """
        updated_sample = insert_size.initialize_sample_attributes(self.sample)
        assert updated_sample.quast.total_reads == 0
        assert updated_sample.quast.insert_mean == []
        assert updated_sample.quast.insert_std == []
        assert updated_sample.quast.read_blocks == []

    @patch(
        'builtins.open', new_callable=mock_open,
        read_data=(
            '# candidate unique pairs for (FF, FR, RF, RR): '
            '(46, 226102, 14, 28)\n'
            'analyzing insert size distribution for orientation FR\n'
            '[M::mem_pestat] mean and std.dev: (487.88, 246.14)\n'
        )
    )
    @patch('cowbat.insert_size.extract_insert_size_distribution')
    def test_parse_reads_stats_file(
        self, mock_extract_insert_size_distribution, mock_open_instance
    ):
        """
        Test the parse_reads_stats_file function.
        """
        mock_extract_insert_size_distribution.return_value = self.sample

        updated_sample = insert_size.parse_reads_stats_file(self.sample)
        mock_open_instance.assert_called_once_with(
            self.sample.quast.reads_stats_file, 'r', encoding='utf-8'
        )
        assert updated_sample.quast.total_reads == 226102
        mock_extract_insert_size_distribution.assert_called_once_with(
            mock_open_instance(), self.sample, 226102
        )

    @patch(
        'builtins.open', new_callable=mock_open,
        read_data='[M::mem_pestat] mean and std.dev: (487.88, 246.14)\n'
    )
    def test_extract_insert_size_distribution(self, mock_open_instance):
        """
        Test the extract_insert_size_distribution function.
        """
        read_stats = mock_open_instance()
        current_reads = 226102
        updated_sample = insert_size.extract_insert_size_distribution(
            read_stats, self.sample, current_reads
        )
        assert updated_sample.quast.insert_mean == [487.88]
        assert updated_sample.quast.insert_std == [246.14]
        assert updated_sample.quast.read_blocks == [226102]

    def test_calculate_weighted_insert_size(self):
        """
        Test the calculate_weighted_insert_size function.
        """
        self.sample.quast.insert_mean = [487.88]
        self.sample.quast.insert_std = [246.14]
        self.sample.quast.read_blocks = [226102]
        self.sample.quast.total_reads = 226102

        updated_metadata = insert_size.calculate_weighted_insert_size(
            self.metadata
        )
        assert self.sample.quast.mean_insert == 487.88
        assert self.sample.quast.std_insert == 246.14
        assert updated_metadata == self.metadata
