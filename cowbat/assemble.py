#!/usr/env/python3

"""
Collection of functions to assemble bacterial genomes
"""

# Standard imports
import logging
from typing import Any, List

# Third-party imports
from cowbat.assembly_evaluation import AssemblyEvaluation
from cowbat.multiqc import multi_qc
from cowbat.prodigal import Prodigal
from cowbat.skesa import Skesa
from olctools.accessoryFunctions.accessoryFunctions import (
    write_metadata_to_file
)

__author__ = 'adamkoziol'


def assemble(
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[Any],
    report_path: str,
    sequence_path: str,
    threads: int
) -> List[Any]:
    """
    Assemble genomes and perform basic quality analyses.

    This function orchestrates the genome assembly process, performs quality
    analyses, and detects open reading frames (ORFs). It updates the metadata
    at each step and writes it to a file.

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
    try:
        # Assemble genomes
        assembly = Skesa(
            log_file=log_file,
            logger=logger,
            metadata=metadata,
            report_path=report_path,
            sequence_path=sequence_path,
            threads=threads
        )
        metadata = assembly.main()

        # Write the metadata to file
        write_metadata_to_file(
            error_logger=error_logger,
            logger=logger,
            metadata=metadata
        )

        # Calculate assembly metrics on raw assemblies
        qual = AssemblyEvaluation(
            log_file=log_file,
            logger=logger,
            metadata=metadata,
            sequence_path=sequence_path,
            threads=threads
        )
        metadata = qual.main()

        # Write the metadata to file
        write_metadata_to_file(
            error_logger=error_logger,
            logger=logger,
            metadata=metadata
        )

        # ORF detection
        prod = Prodigal(
            log_file=log_file,
            logger=logger,
            metadata=metadata
        )
        metadata = prod.main()

        # Write the metadata to file
        write_metadata_to_file(
            error_logger=error_logger,
            logger=logger,
            metadata=metadata
        )

        # Use MultiQC to aggregate reports
        multi_qc(
            log_file=log_file,
            logger=logger,
            report_path=report_path,
            sequence_path=sequence_path
        )

        return metadata

    except (IOError, OSError) as file_error:
        error_logger.error("File operation error: %s", file_error)
        raise
    except RuntimeError as runtime_error:
        error_logger.error(
            "Runtime error during processing: %s", runtime_error
        )
        raise
    except Exception as exc:
        error_logger.error("An unexpected error occurred: %s", exc)
        raise
