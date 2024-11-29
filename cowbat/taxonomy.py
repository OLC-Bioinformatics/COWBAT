#! /usr/env/python3

"""
Collection of methods to perform taxonomy analyses on raw sequence data
and assemblies
"""

# Standard imports
import logging
from typing import List

# Third-party imports
from genemethods.assemblypipeline.mash import run_mash_analyses
from olctools.accessoryFunctions.accessoryFunctions import CustomBox

# Local imports
from cowbat.metaphlan import run_metaphlan_analyses

__author__ = 'adamkoziol'


def taxonomy(
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    report_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Perform taxonomy analyses on raw sequence data and assemblies

    will set sample.general.reference_genus and
    sample.general.closest_refseq_genus to the genus of the closest reference

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[Any]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.
        threads (int): Number of threads to use for processing.

    """
    # Perform Mash analyses
    metadata = run_mash_analyses(
        analysis_type='mash',
        error_logger=error_logger,
        log_file=log_file,
        metadata=metadata,
        reference_file_path=reference_file_path,
        report_path=report_path,
        threads=threads
    )

    # Perform metaphlan analyses on FASTQ files
    metadata = run_metaphlan_analyses(
        error_logger=error_logger,
        file_format='fastq',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        report_path=report_path,
        threads=threads
    )

    # Perform metaphlan analyses on FASTA files
    metadata = run_metaphlan_analyses(
        error_logger=error_logger,
        file_format='fasta',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        report_path=report_path,
        threads=threads
    )

    # Return the updated metadata
    quit()

    return metadata
