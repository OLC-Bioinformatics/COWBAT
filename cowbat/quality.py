#!/usr/bin/env python3

"""
Quality checking functions. Includes FastQC, FASTQ read trimming, error
correction, and contamination discovery with ConFindr.
"""

# Standard imports
import logging
from typing import List

# Third-party imports
from cowbat.contamination_detection import (
    contamination_finder
)
from cowbat.error_correct import error_correction
from cowbat.fastqc import fastqc_threader
from cowbat.quality_trim import trim_quality
from cowbat.repair_gzip import repair_gzip
from cowbat.validate_fastq import validate_fastq
from olctools.accessoryFunctions.accessoryFunctions import (
    write_metadata_to_file
)


def quality(
    *,  # Enforce keyword arguments
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List,
    report_path: str,
    sequence_path: str,
    threads: int
) -> List:
    """
    Run quality checking, trimming, error correction, and contamination
    detection on the samples.

    Args:
        error_logger (logging.Logger): Logger for error messages.
        log_file (str): Path to the log file.
        metadata (List): List of metadata sample objects.
        report_path (str): Path to the report directory.
        sequence_path (str): Path to the sequence directory.
        threads (int): Number of threads to use.

    Returns:
        List: Updated metadata sample objects.
    """
    # Validate that the FASTQ files are in the proper format, and that
    # there are no issues e.g. different numbers of forward and reverse
    # reads, read length longer than quality score length, proper extension
    metadata = validate_fastq(
        error_logger=error_logger,
        metadata=metadata,
        log_file=log_file,
        logger=logger,
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Run FastQC on the unprocessed fastq files
    metadata = fastqc_threader(
        level='Raw',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Perform quality trimming and FastQC on the trimmed files
    metadata = trim_quality(
        log_file=log_file,
        logger=logger,
        metadata=metadata,
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Run FastQC on the trimmed files
    metadata = fastqc_threader(
        level='trimmed',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Perform error correction on the reads
    metadata = error_correction(
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Run FastQC on the processed fastq files
    metadata = fastqc_threader(
        level='trimmed_corrected',
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Fix issue with bbmap gzip
    repair_gzip(
        log_file=log_file,
        logger=logger,
        metadata=metadata
    )

    # Detect contamination in the reads
    metadata = contamination_finder(
        input_path=sequence_path,
        log_file=log_file,
        logger=logger,
        metadata=metadata,
        report_path=report_path,
        threads=threads
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    return metadata
