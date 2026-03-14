#!/usr/env/python3

"""
Collection of functions to performing typing on bacterial sequences
"""

# Standard imports
import logging
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    write_metadata_to_file
)
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.rmlst import rmlst

__author__ = 'adamkoziol'


def agnostic_typing(
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    report_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Perform genus-specific and agnostic typing of bacterial sequences

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[Any]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.
        threads (int): Number of threads to use for the analyses

    Returns:
        List[Any]: Updated metadata after all processing steps.

    Raises:
        IOError: If there is an issue with file operations.
        RuntimeError: If there is an issue with the assembly or quality
            analyses.
    """
    metadata = rmlst(
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
    quit()
    # return metadata
