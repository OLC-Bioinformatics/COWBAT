#!/usr/bin/env python3

"""
FASTQ quality trimming
"""

# Standard imports
from glob import glob
import logging
import os
from subprocess import CalledProcessError
from typing import List

# Third-party imports
from genewrappers.biotools import bbtools
from olctools.accessoryFunctions.accessoryFunctions import write_to_log_file
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def trim_quality(
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Uses bbduk from the bbmap tool suite to quality and adapter trim FASTQ
    files.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info("Trimming FASTQ files")

    # Iterate through samples with FASTQ files
    for sample in metadata:
        logger.debug('Processing sample: %s', sample.name)

        # Only process if .fastq_files is a list
        if not isinstance(sample.general.fastq_files, list):
            logger.debug(
                'Skipping sample %s as fastq_files is not a list', sample.name
            )
            continue

        # Check if the FASTQ files exist
        fastq_files = sorted(sample.general.fastq_files)

        # Define the output directory
        output_dir = sample.general.output_directory

        # Define the names of the trimmed FASTQ files
        clean_forward = os.path.join(
            output_dir, f'{sample.name}_R1_trimmed.fastq.gz'
        )
        clean_reverse = os.path.join(
            output_dir, f'{sample.name}_R2_trimmed.fastq.gz'
        )

        # Set minlength parameter based on read length, default to 50
        try:
            lesser_length = min(
                int(sample.run.Reads.forward_read_length),
                int(sample.run.Reads.reverse_read_length)
            )
        except ValueError:
            lesser_length = int(sample.run.Reads.forward_read_length)
        min_len = 50 if lesser_length >= 50 else lesser_length

        # Initialize variable to store the number of bases to trim
        trim_left = 0

        # Determine if only reverse reads are present
        try:
            if 'R2' in fastq_files[0]:
                if not os.path.isfile(clean_reverse):
                    logger.debug(
                        'Running bbduk_trim for reverse reads of sample: %s',
                        sample.name
                    )
                    out, err, bbduk_call = bbtools.bbduk_trim(
                        forward_in=fastq_files[0],
                        reverse_in=None,
                        forward_out=clean_reverse,
                        trimq=10,
                        minlength=min_len,
                        forcetrimleft=trim_left,
                        returncmd=True
                    )
                else:
                    logger.debug(
                        'Clean reverse file already exists for sample: %s',
                        sample.name
                    )
                    bbduk_call, out, err = '', '', ''
            else:
                if not os.path.isfile(clean_forward):
                    logger.debug(
                        'Running bbduk_trim for forward reads of sample: %s',
                        sample.name
                    )
                    out, err, bbduk_call = bbtools.bbduk_trim(
                        forward_in=fastq_files[0],
                        forward_out=clean_forward,
                        trimq=10,
                        minlength=min_len,
                        forcetrimleft=trim_left,
                        returncmd=True
                    )
                else:
                    logger.debug(
                        'Clean forward file already exists for sample: %s',
                        sample.name
                    )
                    bbduk_call, out, err = '', '', ''
        except (IndexError, CalledProcessError, AttributeError) as exc:
            logger.error(
                'Error during bbduk_trim for sample: %s with error: %s',
                sample.name, str(exc)
            )
            bbduk_call, out, err = '', '', ''

        # Write the command, stdout, and stderr to the log_file
        write_to_log_file(
            out=bbduk_call,
            err=bbduk_call,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )

        # Add the trimmed FASTQ files to a list
        trimmed_fastq_files = sorted(
            glob(os.path.join(output_dir, '*trimmed.fastq.gz'))
        )

        # Populate the metadata if the files exist
        sample.general.trimmed_fastq_files = (
            trimmed_fastq_files if trimmed_fastq_files else []
        )

    logger.info('FASTQ files trimmed')
    return metadata
