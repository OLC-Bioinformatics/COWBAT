#!/usr/bin/env python3

"""
Common methods used by the COWBAT pipeline.
"""

# Standard imports
import logging
import os
from typing import List, Tuple

# Third party imports
from olctools.accessoryFunctions.metadata import CustomBox


def read_checkpoint(checkpoint_file: str) -> str:
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


def write_checkpoint(checkpoint_file: str, step: str) -> None:
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


def tilde_expand(path: str) -> str:
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


def initialize_logging(
    error_log_file: str, log_file: str, logging_level: str
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
    # Create a logger
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, logging_level))

    # Create a formatter without milliseconds
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Check if handlers already exist to avoid duplicate logs
    if not logger.handlers:
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


def initialize_error_log(error_log_file: str) -> logging.Logger:
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

    # Check if error logger handlers already exist
    if not error_logger.handlers:
        # File handler for logging errors to a file
        error_file_handler = logging.FileHandler(error_log_file)
        error_file_handler.setLevel(logging.ERROR)
        error_file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'))
        error_logger.addHandler(error_file_handler)

        # Stream handler for logging errors to the console
        error_console_handler = logging.StreamHandler()
        error_console_handler.setLevel(logging.ERROR)
        error_console_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'))
        error_logger.addHandler(error_console_handler)

    return error_logger


def write_metadata_to_file(
        error_logger: logging.Logger,
        metadata: List[CustomBox]) -> None:
    """
    Write metadata for each sample to its respective JSON file.

    This function iterates over a list of CustomBox metadata objects for
    samples and writes each one to its specified JSON file. If an IOError
    occurs during the file writing process, it logs the error and calls an
    error logger.

    Args:
        error_logger (logging.Logger): Logger instance for logging errors.
        metadata (List[CustomBox]): List of CustomBox metadata objects.

    Example usage:
    >>> metadata = [CustomBox(name='sample1',
        json_file='sample1_metadata.json')]
    >>> write_metadata_to_file(error_logger, metadata)
    """
    for sample in metadata:
        try:
            # Write the metadata to the specified JSON file
            sample.to_file(file_path=sample.json_file)
        except IOError:
            # Log the error traceback
            logging.error("Failed to write metadata to file", exc_info=True)
            # Call the error logger with the error traceback
            error_logger.error(
                "Failed to write metadata to file",
                exc_info=True)
