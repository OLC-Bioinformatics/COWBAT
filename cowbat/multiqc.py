#!/usr/env/python 3

"""
Perform MultiQC analyses on the FASTQ files and assemblies
"""

# Standard imports
import logging
from typing import Tuple

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file,
)


def multi_qc(
    log_file: str,
    logger: logging.Logger,
    report_path: str,
    sequence_path: str
):
    """
    Run MultiQC on the FASTQ files and the assemblies

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        report_path (str): Path to save the report.
        sequence_path (str): Path to the sequence files.
    """
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


def _prepare_command(
    report_path: str,
    sequence_path: str
) -> str:
    """
    Create the MultiQC system call
    """
    command = (
        f'multiqc {sequence_path} --outdir {report_path} --force '
        '--ignore "trimmed*" --ignore-samples "stdin"'
    )

    return command


def _run_command(
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
