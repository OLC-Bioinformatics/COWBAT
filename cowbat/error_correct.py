#!/usr/bin/env python3

"""
FASTQ error correction
"""

# Standard imports
import logging
import os
from subprocess import CalledProcessError
from typing import List

# Third-party imports
from genewrappers.biotools import bbtools
from olctools.accessoryFunctions.accessoryFunctions import write_to_log_file
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def error_correction(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Use tadpole from the bbmap suite of tools to perform error correction
    of the reads.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.
        threads (int): Number of threads to use in the analyses

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Error correcting reads')

    for sample in metadata:
        logger.debug('Processing sample: %s', sample.name)

        # Define the paths for the trimmed corrected FASTQ files
        sample.general.trimmed_corrected_fastq_files = [
            f"{fastq.split('.fastq.gz')[0]}_trimmed_corrected.fastq.gz"
            for fastq in sorted(sample.general.fastq_files)
        ]
        if os.path.isfile(sample.general.trimmed_corrected_fastq_files[0]):
            logger.debug(
                'Trimmed corrected FASTQ file already exists for sample: %s',
                sample.name
            )
            continue

        try:
            logger.debug('Running tadpole for sample: %s', sample.name)
            out, err, cmd = bbtools.tadpole(
                forward_in=sorted(sample.general.trimmed_fastq_files)[0],
                forward_out=sample.general.trimmed_corrected_fastq_files[0],
                returncmd=True,
                mode='correct',
                threads=threads
            )

            # Set the command in the object
            sample.commands.error_correct_cmd = cmd
            logger.debug('Tadpole command: %s', cmd)

            # Write the command, stdout, and stderr to the logfile
            write_to_log_file(
                out=out,
                err=err,
                log_file=log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err,
            )
        except CalledProcessError as exc:
            logger.error(
                'CalledProcessError for sample: %s with error: %s',
                sample.name, str(exc)
            )
            sample.general.trimmed_corrected_fastq_files = (
                sample.general.trimmed_fastq_files
            )
        except AttributeError as exc:
            logger.error(
                'AttributeError for sample: %s with error: %s',
                sample.name, str(exc)
            )
            sample.general.trimmed_corrected_fastq_files = []
        except IndexError as exc:
            logger.error(
                'IndexError for sample: %s with error: %s',
                sample.name, str(exc)
            )
            sample.general.trimmed_corrected_fastq_files = []

    logger.info('Error correction completed')
    return metadata
