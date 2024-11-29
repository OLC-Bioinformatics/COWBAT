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

__author__ = 'adamkoziol'


def typing(
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    report_path: str,
    sequence_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Perform genus-specific and agnostic typing of bacterial sequences

    

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[Any]): List of metadata objects for the samples.
        report_path (str): Path to save the report.
        sequence_path (str): Path to the sequence files.
        threads (int): Number of threads to use for processing.

    Returns:
        List[Any]: Updated metadata after all processing steps.

    Raises:
        IOError: If there is an issue with file operations.
        RuntimeError: If there is an issue with the assembly or quality
            analyses.
    """
    for sample in metadata:
        print(sample.general.reference_genus)
    quit()
