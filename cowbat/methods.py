#!/usr/bin/env python3

"""
Common methods used by the COWBAT pipeline.
"""

# Standard imports
import logging
import os
import subprocess
from typing import (
    List,
    Tuple
)

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    write_metadata_to_file,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.fastq_mover import FastqMover
from cowbat.set_run_metadata import (
    basic,
    determine_and_parse_sample_sheet,
    process_sample
)
from cowbat.teacup_version import __version__


def read_checkpoint(*, checkpoint_file: str) -> str:
    """
    Read the last successful step from the checkpoint file.

    :param checkpoint_file: Path to the checkpoint file
    :return: The last successful step as a string, or an empty string if the
        file does not exist

    Example usage:
    >>> step = read_checkpoint('.checkpoint')
    >>> print(step)  # Output: 'step1' or ''
    """
    # Check if the checkpoint file exists
    if os.path.exists(checkpoint_file):
        # Open the file and read the last successful step
        with open(checkpoint_file, 'r', encoding='utf-8') as f:
            return f.read().strip()
    # Return an empty string if the file does not exist
    return ""


def write_checkpoint(*, checkpoint_file: str, step: str) -> None:
    """
    Write the current step to the checkpoint file.

    :param checkpoint_file: Path to the checkpoint file
    :param step: The current step to be written to the file

    Example usage:
    >>> write_checkpoint('.checkpoint', 'step1')
    """
    # Open the file in write mode and write the current step
    with open(checkpoint_file, 'w', encoding='utf-8') as f:
        f.write(step)


def tilde_expand(*, path: str) -> str:
    """
    Expand the tilde (~) in the supplied path to the full user directory path.

    :param path: The path that may contain a tilde
    :return: The expanded absolute path

    Example usage:
    >>> expanded_path = tilde_expand('~/my_dir')
    >>> print(expanded_path)  # Output: '/home/user/my_dir'
    """
    # Check if the path starts with a tilde (~)
    if path.startswith('~'):
        # Expand the tilde to the full user directory path
        return_path = os.path.abspath(
            os.path.expanduser(
                os.path.join(path)
            )
        )
    else:
        # Return the absolute path if no tilde is present
        return_path = os.path.abspath(
            os.path.join(path)
        )

    return return_path


def initialize_error_log(*, error_log_file: str) -> logging.Logger:
    """
    Initialize the error log file by setting up the file handler.

    Args:
        error_log_file (str): Path to the error log file.

    Returns:
        logging.Logger: The error logger instance.
    """
    # Create a separate logger for errors
    error_logger = logging.getLogger('error_logger')
    error_logger.setLevel(logging.ERROR)

    # Clear existing handlers
    if error_logger.hasHandlers():
        error_logger.handlers.clear()

    # File handler for logging errors to a file
    error_file_handler = logging.FileHandler(error_log_file)
    error_file_handler.setLevel(logging.ERROR)
    error_file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    error_logger.addHandler(error_file_handler)

    # Stream handler for logging errors to the console
    error_console_handler = logging.StreamHandler()
    error_console_handler.setLevel(logging.ERROR)
    error_console_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    error_logger.addHandler(error_console_handler)

    return error_logger


def initialize_logging(
    *, error_log_file: str, log_file: str, logging_level: str
) -> Tuple[logging.Logger, logging.Logger]:
    """
    Initializes logging for the application.

    Sets up logging to both a file and the console with the specified
    logging level. Also sets up a separate error logger.

    Args:
        error_log_file (str): Path to the error log file.
        log_file (str): Path to the main log file.
        logging_level (str): The logging level as a string (e.g., 'DEBUG',
            'INFO').

    Returns:
        Tuple[logging.Logger, logging.Logger]: The main logger instance and
            the error logger instance.
    """
    # Remove previous versions of the files before initializing
    for file in [log_file, error_log_file]:
        if os.path.exists(file):
            os.remove(file)

    # Create a logger
    logger = logging.getLogger('main_logger')
    logger.setLevel(getattr(logging, logging_level))

    # Clear existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create a formatter without milliseconds
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Create a file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(getattr(logging, logging_level))
    file_handler.setFormatter(formatter)

    # Create a stream handler (console)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, logging_level))
    console_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # Initialize the error log
    error_logger = initialize_error_log(error_log_file=error_log_file)

    return logger, error_logger


def sample_metadata(
    *,  # Enforce the use of keyword arguments
    error_logger: logging.Logger,
    logger: logging.Logger,
    metadata: List[CustomBox],
    sequence_path: str
):
    """
    Helper method to prep metadata from sample sheet (if available)

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata objects for the samples.
        sequence_path (str): Path to the sequence files.
    """
    # Define the sample sheet
    sample_sheet = os.path.join(
        sequence_path,
        'SampleSheet.csv'
    )

    # Process the samples appropriate depending on whether a sample sheet
    # was provided
    if os.path.isfile(sample_sheet):

        # Extract the necessary information from the sample sheet
        data = determine_and_parse_sample_sheet(
            file_path=sample_sheet,
            logger=error_logger
        )

        # Create the metadata for each sample
        metadata = process_sample(
            commit=__version__,
            data=data,
            logger=error_logger,
            path=sequence_path
        )
    else:
        # If the sample sheet is missing, perform basic assembly
        metadata = basic(
            commit=__version__,
            logger=error_logger,
            sequence_path=sequence_path
        )
    # Move/link the FASTQ files to strain-specific working directories
    FastqMover(
        metadata=metadata,
        logger=logger,
        path=sequence_path
    )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    return metadata


def check_programs(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    programs: List[str] = None
) -> List[str]:
    """
    Check if the required command line programs are installed.

    Args:
        logger (logging.Logger): Logger for recording information.
        programs (List[str]): List of program names to check. Defaults to None.

    Returns:
        List[str]: List of missing programs.
    """
    if not programs:
        logger.error("No programs provided to check.")
        return []

    missing_programs = []

    for program in programs:
        try:
            # Attempt to run the program with a simple command to check if it
            # is installed
            subprocess.run(
                [program, '--version'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
            logger.debug("%s is installed.", program)
        except subprocess.CalledProcessError:
            try:
                # Attempt to run the program with a different command
                subprocess.run(
                    [program, '-h'],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True
                )
                logger.debug("%s is installed.", program)
            except subprocess.CalledProcessError:
                logger.warning("%s is missing.", program)
                missing_programs.append(program)
        except FileNotFoundError:
            logger.warning("%s is missing.", program)
            missing_programs.append(program)

    if missing_programs:
        logger.error("The following programs are missing:")
        for program in missing_programs:
            logger.error(program)
        return missing_programs

    logger.info("All required programs are installed.")
    return []


def write_to_log_files(
    *,  # Enforce keyword arguments
    command: str,
    err: str,
    log_file: str,
    logger: logging.Logger,
    out: str,
    program: str,
    sample: CustomBox,
    log_output: bool = False
) -> None:
    """
    Write the command and outputs to the log files

    Args:
        command (str): The command.
        err (str): The error output from the command.
        file_format (str): The format of the sequence files.
        log_file (str): The name and path of the log file.
        logger (logging.Logger): The logger object.
        out (str): The output from the command.
        program (str): The name of the program for which the command was run.
        sample (CustomBox): The sample object
        log_output (bool): Whether to log the output to the sample log file.
    """
    # Log the command
    if log_output:
        logger.debug(
            'Writing %s command to log for sample: %s',
            program, sample.name
        )

    # Write the command to the log files
    write_to_log_file(
        out=command,
        err=command,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err,
    )

    # Log the output from the command
    if log_output:
        logger.debug(
            'Writing %s output to log for sample: %s',
            program, sample.name
        )

    # Write the outputs to the log files
    write_to_log_file(
        out=out,
        err=err,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err,
    )
