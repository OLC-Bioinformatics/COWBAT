# tests/test_quality_report.py

"""
Unit tests for quality report functions in the Teacup COWBAT pipeline
"""

# Standard imports
import os
from datetime import datetime
from unittest.mock import mock_open, patch

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Import functions directly from cowbat.quality_report
from cowbat.quality_report import (
    _determine_download_date,
    extract_attribute,
    extract_n50,
    generate_sample_data,
    write_quality_report
)

# Import version from cowbat.teacup_version
from cowbat.teacup_version import __version__


class TestQualityReport:
    """
    Unit tests for Teacup COWBAT quality report functions
    """

    @classmethod
    def setup_class(cls):
        """
        Setup method to initialize any state before tests run.
        """
        cls.commit = __version__
        cls.reference_file_path = 'tests/testdata/reference'
        cls.report_path = 'tests/testdata/report'
        cls.metadata = [
            CustomBox({
                'name': 'sample1',
                'run': CustomBox(Data=CustomBox(Sample_Plate='Plate1')),
                'general': CustomBox(closest_refseq_genus='Genus1'),
                'confindr': CustomBox(num_contaminated_snvs=5),
                'quast': CustomBox(
                    N50=5000,
                    num_contigs=10,
                    Total_length=50000,
                    mean_insert=300,
                    std_insert=30,
                    GC=50
                ),
                'qualimap': CustomBox(
                    MeanCoveragedata=30,
                    StdCoveragedata=5
                ),
                'mash': CustomBox(
                    closest_refseq='RefSeq1',
                    num_matches=100
                ),
                'prodigal': CustomBox(
                    predicted_genes_total=4000,
                    predicted_genes_over_3000bp=10,
                    predicted_genes_over_1000bp=50,
                    predicted_genes_over_500bp=100,
                    predicted_genes_under_500bp=3840
                )
            })
        ]

    def test_extract_attribute(self):
        """
        Test the extract_attribute function.
        """
        sample = self.metadata[0]
        assert extract_attribute(sample, 'name') == 'sample1\t'
        assert extract_attribute(
            sample.run.Data, 'Sample_Plate'
        ) == 'Plate1\t'
        assert extract_attribute(
            sample.general, 'closest_refseq_genus'
        ) == 'Genus1\t'
        assert extract_attribute(
            sample.confindr, 'num_contaminated_snvs', number=True
        ) == '5\t'
        assert extract_attribute(
            sample.quast, 'N50', number=True
        ) == '5000\t'
        assert extract_attribute(
            sample.quast, 'num_contigs', number=True
        ) == '10\t'
        assert extract_attribute(
            sample.quast, 'Total_length', number=True
        ) == '50000\t'
        assert extract_attribute(
            sample.quast, 'mean_insert', number=True
        ) == '300\t'
        assert extract_attribute(
            sample.quast, 'std_insert', number=True
        ) == '30\t'
        assert extract_attribute(
            sample.qualimap, 'MeanCoveragedata', number=True
        ) == '30\t'
        assert extract_attribute(
            sample.qualimap, 'StdCoveragedata', number=True
        ) == '5\t'
        assert extract_attribute(
            sample.quast, 'GC', number=True
        ) == '50\t'
        assert extract_attribute(
            sample.mash, 'closest_refseq'
        ) == 'RefSeq1\t'
        assert extract_attribute(
            sample.mash, 'num_matches', number=True
        ) == '100\t'
        assert extract_attribute(
            sample.prodigal, 'predicted_genes_total', number=True
        ) == '4000\t'
        assert extract_attribute(
            sample.prodigal, 'predicted_genes_over_3000bp', number=True
        ) == '10\t'
        assert extract_attribute(
            sample.prodigal, 'predicted_genes_over_1000bp', number=True
        ) == '50\t'
        assert extract_attribute(
            sample.prodigal, 'predicted_genes_over_500bp', number=True
        ) == '100\t'
        assert extract_attribute(
            sample.prodigal, 'predicted_genes_under_500bp', number=True
        ) == '3840\t'

    def test_extract_n50(self):
        """
        Test the extract_n50 function.
        """
        sample = self.metadata[0]
        assert extract_n50(sample) == '5000\t'
        sample.quast.N50 = 0
        assert extract_n50(sample) == '0\t'

    @patch('builtins.open', new_callable=mock_open, read_data='2023-01-01')
    def test_determine_download_date(self, mock_open_instance):
        """
        Test the _determine_download_date function.
        """
        download_date = _determine_download_date(
            self.reference_file_path
        )
        mock_open_instance.assert_called_once_with(
            os.path.join(
                self.reference_file_path, 'download_date'
            ),
            'r', encoding='utf-8'
        )
        assert download_date == '2023-01-01'

    @patch('cowbat.quality_report._determine_download_date')
    def test_generate_sample_data(self, mock_determine_download_date):
        """
        Test the generate_sample_data function.
        """
        mock_determine_download_date.return_value = '2023-01-01'
        sample = self.metadata[0]
        sample.quast.N50 = 5000  # Ensure N50 is set correctly
        sample_data = generate_sample_data(
            sample=sample,
            commit=self.commit,
            reference_file_path=self.reference_file_path
        )
        expected_data = (
            'sample1\tPlate1\tGenus1\t5\t5000\t10\t50000\t300\t30\t30\t5\t50\t'
            f'RefSeq1\t100\t4000\t10\t50\t100\t3840\t{__version__}\t'
            f'{datetime.now().strftime("%Y-%m-%d")}\t'
            'reference\t2023-01-01\n'
        )
        assert sample_data == expected_data

    @patch('builtins.open', new_callable=mock_open)
    @patch('cowbat.quality_report.generate_sample_data')
    def test_write_quality_report(
        self, mock_generate_sample_data, mock_open_instance
    ):
        """
        Test the write_quality_report function.
        """
        mock_generate_sample_data.return_value = (
            'sample1\tPlate1\tGenus1\t5\t5000\t10\t50000\t300\t30\t30\t5\t50\t'
            f'RefSeq1\t100\t4000\t10\t50\t100\t3840\t{__version__}\t'
            f'{datetime.now().strftime("%Y-%m-%d")}\t'
            'reference\t2023-01-01\n'
        )
        write_quality_report(
            metadata=self.metadata,
            commit=self.commit,
            reference_file_path=self.reference_file_path,
            report_path=self.report_path
        )
        mock_open_instance.assert_called_once_with(
            os.path.join(
                self.report_path, 'quality_report.tsv'
            ),
            'w', encoding='utf-8'
        )
        mock_open_instance().write.assert_any_call(
            'SeqID\tSample_Name\tGenus\tConfindr_Contaminated_SNVs\tN50\t'
            'Num_Contigs\tTotal_Length\tMean_Insert_Size\tInsert_Size_STD\t'
            'Average_Coverage_Depth\tCoverage_Depth_STD\tPercent_GC\t'
            'MASH_Reference_Genome\tMASH_Num_Matching_Hashes\t'
            'Total_Predicted_Genes\tPredicted_Genes_Over_3000bp\t'
            'Predicted_Genes_Over_1000bp\tPredicted_Genes_Over_500bp\t'
            'Predicted_Genes_Under_500bp\tAssembly_Date\tPipeline_Version\t'
            'Database\tDatabase_Download_Date\n'
        )
        mock_open_instance().write.assert_any_call(
            'sample1\tPlate1\tGenus1\t5\t5000\t10\t50000\t300\t30\t30\t5\t50\t'
            f'RefSeq1\t100\t4000\t10\t50\t100\t3840\t{__version__}\t'
            f'{datetime.now().strftime("%Y-%m-%d")}\t'
            'reference\t2023-01-01\n'
        )
