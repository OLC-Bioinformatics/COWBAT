#!/usr/env/python 3

"""
Perform MultiQC analyses on the FASTQ files and assemblies
"""

# Standard imports
import logging
from typing import (
    List, Tuple
)

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file,
)
from olctools.accessoryFunctions.metadata import CustomBox


def multi_qc(
    *,  # Enforce the use of keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    report_path: str,
    sequence_path: str
):
    """
    Run MultiQC on the outputs of the pipeline

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata objects for the samples.
        report_path (str): Path to save the report.
        sequence_path (str): Path to the sequence files.
    """
    # Clean the log file to remove invalid literals for int() with base 10
    _clean_log_files(
        metadata=metadata,
    )

    # Prepare the MultiQC command
    command = _prepare_command(
        report_path=report_path,
        sequence_path=sequence_path
    )

    # Run the command
    out, err = _run_command(
        command=command,
        logger=logger
    )

    # Write the outputs to the log file
    _write_to_logfile(
        err=err,
        log_file=log_file,
        out=out
    )

    # Return the metadata
    return metadata


def _clean_log_files(
    *,  # Enforce the use of keyword arguments
    metadata: List[CustomBox]
) -> None:
    """
    Clean the log files to remove invalid literals for int() with base 10

    Args:
        metadata (List[CustomBox]): List of metadata objects for the samples.
    """
    for sample in metadata:
        # Read the log file
        with open(sample.general.log_out, 'r', encoding='utf-8') as log_file:
            log_lines = log_file.readlines()

            # Create a list to store the new log lines
            new_log_lines = []

            # Look for the invalid literals
            for line in log_lines:
                # Search for specific strings
                illegals = ['N50', 'N90']

                # Check if the line contains any of the invalid literals
                if any(illegal in line for illegal in illegals):
                    # Remove the comma from the size e.g. 557,735
                    line = line.replace(',', '')
                new_log_lines.append(line)

        # Write the new log lines to the log file
        with open(sample.general.log_out, 'w', encoding='utf-8') as log_file:
            log_file.writelines(new_log_lines)


def _prepare_command(
    *,  # Enforce the use of keyword arguments
    report_path: str,
    sequence_path: str
) -> str:
    """
    Create the MultiQC system call
    """
    command = (
        f'multiqc {sequence_path} --outdir {report_path} --force '
        '--ignore "trimmed*" --ignore-samples "stdin" '
        '--ignore "*_log_out*" --ignore "*_log_err*" '
    )
    print(command)
    return command


def _run_command(
    *,  # Enforce the use of keyword arguments
    command: str,
    logger: logging.Logger
) -> Tuple[str, str]:
    """
    Run a subprocess command and return the output and error.

    Args:
        command (str): The command to run.
        logger (logging.Logger): Logger for recording information.

    Returns:
        tuple: Standard output and error strings.
    """
    logger.debug('Running MultiQC command: %s', command)
    out, err = run_subprocess(command=command)
    return out, err


def _write_to_logfile(
    *,  # Enforce the use of keyword arguments
    err: str,
    log_file: str,
    out: str
) -> None:
    """
    Write the stdout and stderr to the log file
    """
    write_to_log_file(
        out=out,
        err=err,
        log_file=log_file
    )
