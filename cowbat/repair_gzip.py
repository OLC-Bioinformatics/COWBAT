#!/usr/bin/env python3

"""
Occasionally, downstream programs have issues with the gzip performed by
BBduk. Decompress, and recompress to ensure compatibility.
"""

# Standard imports
import logging
import os
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)

__author__ = 'adamkoziol'


def repair_gzip(
    log_file: str,
    logger: logging.Logger,
    metadata: List
) -> None:
    """
    Decompress and recompress the FASTQ file to fix issues linked to
    BBDuk compression.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List): List of metadata sample objects.
    """
    logger.info('Fixing BBMap gunzip issue')

    for sample in metadata:
        logger.debug('Processing sample: %s', sample.name)
        try:
            fastq_files = sample.general.trimmed_corrected_fastq_files
        except AttributeError:
            fastq_files = None

        # Guard statement to check if fastq_files is a list
        if not isinstance(fastq_files, list):
            logger.warning('No FASTQ files found for sample: %s', sample.name)
            continue

        output_directory = sample.general.output_directory
        tmp_forward = os.path.join(
            output_directory, f'{sample.name}_tmp_R1.fastq'
        )
        tmp_reverse = os.path.join(
            output_directory, f'{sample.name}_tmp_R2.fastq'
        )

        # Construct system calls for decompressing and recompressing FASTQ
        # files
        decompress_recompress_forward = (
            f'gunzip -c {fastq_files[0]} > {tmp_forward} && '
            f'gzip -c {tmp_forward} > {fastq_files[0]} && rm {tmp_forward}'
        )
        decompress_recompress_reverse = (
            f'gunzip -c {fastq_files[1]} > {tmp_reverse} && '
            f'gzip -c {tmp_reverse} > {fastq_files[1]} && rm {tmp_reverse}'
        )

        logger.debug(
            'Running system call for forward read: %s',
            decompress_recompress_forward)

        # Run the system call for the forward read
        out, err = run_subprocess(decompress_recompress_forward)

        write_to_log_file(
            out=decompress_recompress_forward,
            err=decompress_recompress_forward,
            log_file=log_file
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file
        )

        logger.debug(
            'Running system call for reverse read: %s',
            decompress_recompress_reverse)

        # Run the system call for the reverse read
        out, err = run_subprocess(decompress_recompress_reverse)
        write_to_log_file(
            out=decompress_recompress_reverse,
            err=decompress_recompress_reverse,
            log_file=log_file
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file
        )

    logger.info('Fixing BBMap gunzip issue complete!')
