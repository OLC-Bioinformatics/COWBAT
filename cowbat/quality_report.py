#!/usr/env/python3

"""
Create a preliminary report for COWBAT with basic quality and Mash outputs
"""

# Standard imports
import logging
import os
from datetime import datetime
from typing import Any, List

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.multiqc import multi_qc

# Constants
QUALITY_HEADERS = [
    'SeqID',
    'Sample_Name',
    'Genus',
    'Confindr_Contaminated_SNVs',
    'N50',
    'Num_Contigs',
    'Total_Length',
    'Mean_Insert_Size',
    'Insert_Size_STD',
    'Average_Coverage_Depth',
    'Coverage_Depth_STD',
    'Percent_GC',
    'MASH_Reference_Genome',
    'MASH_Num_Matching_Hashes',
    'Total_Predicted_Genes',
    'Predicted_Genes_Over_3000bp',
    'Predicted_Genes_Over_1000bp',
    'Predicted_Genes_Over_500bp',
    'Predicted_Genes_Under_500bp',
    'Assembly_Date',
    'Pipeline_Version',
    'Database',
    'Database_Download_Date'
]


def extract_attribute(
    *,  # All parameters below must be passed as keyword arguments
    attribute: str,
    logger: logging.Logger,
    sample: Any,
    number: bool = False
) -> str:
    """
    Extract the value of the specified attribute from the CustomBox object.

    Args:
        attribute (str): The name of the attribute to extract.
        logger (logging.Logger): Logger object.
        number (bool): Whether to format the attribute as a number.
        sample (Any): The CustomBox object to extract the attribute from.

    Returns:
        str: The value of the specified attribute as a string.

    Raises:
        AttributeError: If the attribute does not exist.
    """
    try:
        value = getattr(sample, attribute)
        if number:
            return f"{value}\t" if value is not None else "0\t"
        return f"{value}\t"
    except AttributeError:
        logger.error(
            "Attribute '%s' not found in sample: %s", attribute, sample
        )
        return "NA\t"


def extract_n50(*, logger: logging.Logger, sample: Any) -> str:
    """
    Extract the N50 attribute with special handling for zero values.

    Args:
        logger (logging.Logger): Logger object.
        sample (Any): The CustomBox object to extract the attribute from.

    Returns:
        str: The value of the N50 attribute as a string.
    """
    n50 = extract_attribute(
        attribute='N50',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    return n50 if n50 != '0\t' else '0\t'


def _determine_download_date(*, reference_file_path: str):
    """
    Extract the database download date from file

    Args:
        reference_file_path (str): The absolute path to the folder containing
            the database
    """
    # Set the name of the file containing the download date
    date_file = os.path.join(reference_file_path, 'download_date')

    try:
        with open(date_file, 'r', encoding='utf-8') as download_date:
            date = download_date.readline().rstrip()
    except FileNotFoundError:
        date = 'ND'

    return date


def generate_sample_data(
    *,
    commit: str,
    logger: logging.Logger,
    reference_file_path: str,
    sample: CustomBox,
) -> str:
    """
    Generate the data string for a single sample.

    Args:
        commit (str): The pipeline version.
        logger (logging.Logger): Logger object.
        sample (CustomBox): The sample object containing metadata.
        reference_file_path (str): The path to the reference file.

    Returns:
        str: The data string for the sample.
    """
    data = ""
    data += extract_attribute(
        attribute='name',
        logger=logger,
        sample=sample
    )
    data += extract_attribute(
        attribute='Sample_Plate',
        logger=logger,
        sample=sample.run.Data
    )
    data += extract_attribute(
        attribute='closest_refseq_genus',
        logger=logger,
        sample=sample.general
    )
    data += extract_attribute(
        attribute='num_contaminated_snvs',
        logger=logger,
        sample=sample.confindr
    )
    data += extract_n50(
        logger=logger,
        sample=sample
    )
    data += extract_attribute(
        attribute='num_contigs',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    data += extract_attribute(
        attribute='Total_length',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    data += extract_attribute(
        attribute='mean_insert',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    data += extract_attribute(
        attribute='std_insert',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    data += extract_attribute(
        attribute='MeanCoveragedata',
        logger=logger,
        number=True,
        sample=sample.qualimap,
    )
    data += extract_attribute(
        attribute='StdCoveragedata',
        logger=logger,
        number=True,
        sample=sample.qualimap
    )
    data += extract_attribute(
        attribute='GC',
        logger=logger,
        number=True,
        sample=sample.quast
    )
    data += extract_attribute(
        attribute='closest_refseq',
        logger=logger,
        sample=sample.mash
    )
    data += extract_attribute(
        attribute='num_matches',
        logger=logger,
        sample=sample.mash
    )
    data += extract_attribute(
        attribute='predicted_genes_total',
        logger=logger,
        number=True,
        sample=sample.prodigal
    )
    data += extract_attribute(
        attribute='predicted_genes_over_3000bp',
        logger=logger,
        number=True,
        sample=sample.prodigal
    )
    data += extract_attribute(
        attribute='predicted_genes_over_1000bp',
        logger=logger,
        number=True,
        sample=sample.prodigal,
    )
    data += extract_attribute(
        attribute='predicted_genes_over_500bp',
        logger=logger,
        number=True,
        sample=sample.prodigal,
    )
    data += extract_attribute(
        attribute='predicted_genes_under_500bp',
        logger=logger,
        number=True,
        sample=sample.prodigal,
    )
    data += f"{commit}\t"
    data += f"{datetime.now().strftime('%Y-%m-%d')}\t"
    data += f"{os.path.split(reference_file_path)[-1]}\t"
    data += _determine_download_date(
        reference_file_path=reference_file_path
    )
    data += '\n'
    return data


def write_quality_report(
    *,  # Enforce keyword arguments
    commit: str,
    log_file: str,
    logger: logging.Logger,
    metadata: List[Any],
    reference_file_path: str,
    report_path: str,
    sequence_path: str
) -> None:
    """
    Create a sample quality summary report in TSV format.

    This function generates a TSV report with basic quality and Mash outputs
    for each sample in the metadata.

    Args:
        commit (str): The pipeline version.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List[Any]): A list of sample objects containing metadata.
        reference_file_path (str): The path to the reference file.
        report_path (str): The path to save the report.
        sequence_path (str): The path to the sequence files.

    Returns:
        None
    """
    # Use MultiQC to aggregate reports
    multi_qc(
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        report_path=report_path,
        sequence_path=sequence_path
    )

    logger.info('Creating sample quality summary report')
    header = f"{'\t'.join(QUALITY_HEADERS)}\n"

    # Create a string to store all the results
    data = ""
    for sample in metadata:
        data += generate_sample_data(
            commit=commit,
            logger=logger,
            reference_file_path=reference_file_path,
            sample=sample,
        )

    # Replace any NA values with ND
    clean_data = data.replace('NA', 'ND')

    # Set the name of the report
    report = os.path.join(report_path, 'quality_report.tsv')

    with open(report, 'w', encoding='utf-8') as quality_report:
        quality_report.write(header)
        quality_report.write(clean_data)

    # Return the metadata
    return metadata
