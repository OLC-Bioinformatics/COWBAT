#! /usr/env/python3

"""
Collection of methods to perform taxonomy analyses on raw sequence data
and assemblies
"""

# Standard imports
import logging
import os
from typing import List

# Third-party imports
from genemethods.assemblypipeline.mash import run_mash_analyses
from olctools.accessoryFunctions.accessoryFunctions import (
    CustomBox,
    write_metadata_to_file
)

# Local imports
from cowbat.metaphlan import run_metaphlan_analyses
from cowbat.gtdbtk import (
    gtdbtk,
    parse_gtbdtk_output
)

__author__ = 'adamkoziol'


def taxonomy(
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    report_path: str,
    sequence_path: str,
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
        sequence_path (str): Path to the sequence files.
        threads (int): Number of threads to use for processing.

    """
    # Perform Mash analyses
    metadata = run_mash_analyses(
        analysis_type='mash',
        error_logger=error_logger,
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        reference_file_path=reference_file_path,
        report_path=report_path,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Perform metaphlan analyses on FASTQ files
    metadata = run_metaphlan_analyses(
        error_logger=error_logger,
        file_format='fastq',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        report_path=report_path,
        sequence_path=sequence_path,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Run GTDB-Tk analyses
    report = gtdbtk(
        assembly_path=os.path.join(sequence_path, 'BestAssemblies'),
        log_file=log_file,
        logger=logger,
        reference_file_path=reference_file_path,
        report_path=report_path,
        threads=threads
    )

    # Parse the GTDB-Tk output
    metadata = parse_gtbdtk_output(
        metadata=metadata,
        report=report
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Return the updated metadata
    return metadata
