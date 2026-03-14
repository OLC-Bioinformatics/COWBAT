#!/usr/bin/env python3

"""
FASTQ validation, reformatting, and repair functions
"""

# Standard imports
import logging
import os
from subprocess import CalledProcessError
import traceback
from typing import List

# Third-party imports
from genewrappers.biotools import bbtools
from olctools.accessoryFunctions.accessoryFunctions import (
    write_to_log_file,
    write_metadata_to_file
)
from olctools.accessoryFunctions.metadata import CustomBox


def validate_fastq(
    *,  # Enforce keyword arguments
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Runs reformat.sh on the FASTQ files. If a CalledProcessError arises, do
    not proceed with the assembly of these files.

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Validating FASTQ files')

    # Iterate over each sample in the metadata
    for sample in metadata:
        logger.debug('Validating sample: %s', sample.name)

        # Check if the FASTQ file size is valid
        if not validate_file_size(
            logger=logger,
            sample=sample
        ):
            # Log an error if the file size is too small
            logger.error('File size too small for sample: %s', sample.name)
            error(
                logger=logger,
                message='files_too_small',
                sample=sample
            )
            continue

        try:
            # Validate the reads using reformat.sh
            logger.debug('Running reformat.sh for sample: %s', sample.name)
            validate_reads(
                log_file=log_file,
                logger=logger,
                sample=sample
            )
            logger.debug('Validation successful for sample: %s', sample.name)
        except CalledProcessError as exc:
            # Handle any validation errors by attempting to repair the reads
            logger.error(
                'Validation failed for sample: %s with error: %s',
                sample.name, str(exc)
            )
            handle_validation_error(
                log_file=log_file,
                logger=logger,
                sample=sample
            )

    # Write the updated metadata to a file
    logger.debug('Writing updated metadata to file')
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata,
    )

    logger.info('FASTQ validation completed')
    return metadata


def validate_file_size(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> bool:
    """
    Validate the size of the FASTQ file.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        bool: True if the file size is valid, False otherwise.
    """
    logger.debug('Validating file size for sample: %s', sample.name)

    # Get the size of the first FASTQ file
    fastq_file = sample.general.fastq_files[0]
    size = os.path.getsize(fastq_file)
    logger.debug('File size for %s: %d bytes', fastq_file, size)

    # Check if the file size is greater than or equal to 1,000,000 bytes
    is_valid = size >= 1000000
    if not is_valid:
        logger.warning('File size too small for sample: %s', sample.name)

    return is_valid


def validate_reads(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Validate the reads using reformat.sh.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        sample (CustomBox): Metadata sample object.

    Raises:
        CalledProcessError: If reformat.sh fails.
    """
    logger.debug('Validating reads for sample: %s', sample.name)

    try:
        # Run reformat.sh to validate the reads
        logger.debug('Running reformat.sh for sample: %s', sample.name)
        out, err, _ = bbtools.validate_reads(
            forward_in=sample.general.fastq_files[0],
            returncmd=True
        )
        logger.debug('reformat.sh completed for sample: %s', sample.name)

        # Write the output and error messages to the log file
        logger.debug(
            'Writing reformat.sh output to log file for sample: %s',
            sample.name)
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except CalledProcessError as exc:
        logger.error(
                'reformat.sh failed for sample: %s with error: %s',
                sample.name, str(exc)
            )
        logger.error(traceback.format_exc())
        raise


def handle_validation_error(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Handle validation errors by attempting to repair the reads.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        sample (CustomBox): Metadata sample object.
    """
    # Log a warning about detected errors in the FASTQ files
    logger.warning(
        'Errors detected in FASTQ files for sample %s. '
        'Please check the following files for details %s %s %s. '
        'The pipeline will use reformat.sh to attempt to repair issues',
        sample.name, log_file, sample.general.log_out,
        sample.general.log_err
    )

    # Get the paths for reformatted and repaired files
    reformatted_forward, reformatted_reverse, repair_forward, \
        repair_reverse = get_reformatted_and_repair_paths(
            logger=logger,
            sample=sample
        )
    logger.debug(
        'Reformatted and repair paths for sample %s: %s, %s, %s, %s',
        sample.name, reformatted_forward, reformatted_reverse,
        repair_forward, repair_reverse
    )

    # If the reformatted forward file already exists, return early
    if os.path.isfile(reformatted_forward):
        logger.debug(
            'Reformatted forward file already exists for sample %s: %s',
            sample.name, reformatted_forward
        )
        return

    try:
        # Attempt to reformat the reads
        logger.debug(
            'Attempting to reformat reads for sample %s',
            sample.name)
        reformat_reads(
            log_file=log_file,
            logger=logger,
            sample=sample,
            reformatted_forward=reformatted_forward
        )
        logger.debug('Reformatting completed for sample %s', sample.name)

        # Attempt to repair the reads
        logger.debug('Attempting to repair reads for sample %s', sample.name)
        repair_reads(
            log_file=log_file,
            logger=logger,
            sample=sample,
            reformatted_forward=reformatted_forward,
            reformatted_reverse=reformatted_reverse,
            repair_forward=repair_forward,
            repair_reverse=repair_reverse
        )
        logger.debug('Repairing completed for sample %s', sample.name)
    except CalledProcessError as exc:
        # Handle any errors that occur during the repair process
        logger.error(
            'Repair process failed for sample %s with error: %s',
            sample.name, str(exc)
        )
        logger.error(traceback.format_exc())
        handle_repair_error(
            log_file=log_file,
            logger=logger,
            sample=sample,
            reformatted_forward=reformatted_forward,
            reformatted_reverse=reformatted_reverse,
            repair_forward=repair_forward,
            repair_reverse=repair_reverse
        )


def get_reformatted_and_repair_paths(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> tuple:
    """
    Get the paths for reformatted and repaired files.

    Args:
        logger (logging.Logger): Logger object.
        sample (CustomBox): Metadata sample object.

    Returns:
        tuple: Paths for reformatted and repaired files.
    """
    logger.debug(
        'Getting reformatted and repair paths for sample: %s',
        sample.name)

    # Define the path for the reformatted forward file
    reformatted_forward = os.path.join(
        sample.general.output_directory,
        f'{sample.name}_reformatted_R1.fastq.gz'
    )
    logger.debug('Reformatted forward path: %s', reformatted_forward)

    # Define the path for the repaired forward file
    repair_forward = os.path.join(
        sample.general.output_directory,
        f'{sample.name}_repaired_R1.fastq.gz'
    )
    logger.debug('Repair forward path: %s', repair_forward)

    # Check if there are two FASTQ files
    if len(sample.general.fastq_files) == 2:
        # Define the paths for the reformatted and repaired reverse files
        reformatted_reverse = os.path.join(
            sample.general.output_directory,
            f'{sample.name}_reformatted_R2.fastq.gz'
        )
        repair_reverse = os.path.join(
            sample.general.output_directory,
            f'{sample.name}_repaired_R2.fastq.gz'
        )
        logger.debug('Reformatted reverse path: %s', reformatted_reverse)
        logger.debug('Repair reverse path: %s', repair_reverse)
    else:
        # If there is only one FASTQ file, set the reverse paths to empty
        # strings
        reformatted_reverse = ''
        repair_reverse = ''
        logger.debug('Single FASTQ file detected, no reverse paths needed.')

    # Return the paths for the reformatted and repaired files
    return reformatted_forward, reformatted_reverse, repair_forward, \
        repair_reverse


def reformat_reads(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    reformatted_forward: str,
    sample: CustomBox,
) -> None:
    """
    Reformat the reads using reformat.sh.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        reformatted_forward (str): Path to the reformatted forward file.
        sample (CustomBox): Metadata sample object.

    Raises:
        CalledProcessError: If reformat.sh fails.
    """
    logger.debug('Reformatting reads for sample: %s', sample.name)

    try:
        # Run reformat.sh to reformat the reads
        logger.debug(
            'Running reformat.sh for sample: %s, output: %s',
            sample.name, reformatted_forward
        )
        out, err, _ = bbtools.reformat_reads(
            forward_in=sample.general.fastq_files[0],
            forward_out=reformatted_forward,
            returncmd=True
        )
        logger.debug('reformat.sh completed for sample: %s', sample.name)

        # Write the output and error messages to the log file
        logger.debug(
            'Writing reformat.sh output to log file for sample: %s',
            sample.name
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except CalledProcessError as exc:
        logger.error(
            'reformat.sh failed for sample: %s with error: %s',
            sample.name, str(exc)
        )
        logger.error(traceback.format_exc())
        raise


def repair_reads(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    reformatted_forward: str,
    reformatted_reverse: str,
    repair_forward: str,
    repair_reverse: str,
    sample: CustomBox
) -> None:
    """
    Repair the reads using repair.sh.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        reformatted_forward (str): Path to the reformatted forward file.
        reformatted_reverse (str): Path to the reformatted reverse file.
        repair_forward (str): Path to the repaired forward file.
        repair_reverse (str): Path to the repaired reverse file.
        sample (CustomBox): Metadata sample object.

    Raises:
        CalledProcessError: If repair.sh fails.
    """
    logger.debug('Repairing reads for sample: %s', sample.name)

    try:
        # Check if there is a reformatted reverse file
        if reformatted_reverse:
            logger.debug(
                'Running repair.sh for sample: %s, forward: %s, reverse: %s',
                sample.name, reformatted_forward, reformatted_reverse
            )
            # Run repair.sh to repair the reads
            out, err, _ = bbtools.repair_reads(
                forward_in=reformatted_forward,
                reverse_in=reformatted_reverse,
                forward_out=repair_forward,
                reverse_out=repair_reverse,
                returncmd=True
            )
            logger.debug('repair.sh completed for sample: %s', sample.name)

            # Write the output and error messages to the log file
            logger.debug(
                'Writing repair.sh output to log file for sample: %s',
                sample.name
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

        # Check if the reformatted forward file exists
        if os.path.isfile(reformatted_forward):
            logger.debug(
                'Updating fastq_files attribute for sample: %s', sample.name
            )
            # Update the fastq_files attribute to point to the repaired files
            sample.general.fastq_files = (
                [repair_forward, repair_reverse]
                if repair_reverse else [reformatted_forward]
            )
    except CalledProcessError as exc:
        logger.error(
            'repair.sh failed for sample: %s with error: %s',
            sample.name, str(exc)
        )
        logger.error(traceback.format_exc())
        raise


def handle_repair_error(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    reformatted_forward: str,
    reformatted_reverse: str,
    repair_forward: str,
    repair_reverse: str,
    sample: CustomBox
) -> None:
    """
    Handle errors during the repair process.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        reformatted_forward (str): Path to the reformatted forward file.
        reformatted_reverse (str): Path to the reformatted reverse file.
        repair_forward (str): Path to the repaired forward file.
        repair_reverse (str): Path to the repaired reverse file.
        sample (CustomBox): Metadata sample object.
    """
    logger.debug('Handling repair error for sample: %s', sample.name)

    # Check if both reformatted forward and reverse files exist
    if os.path.isfile(reformatted_forward) and os.path.isfile(
            reformatted_reverse):
        try:
            logger.debug(
                'Running repair.sh for sample: %s, forward: %s, reverse: %s',
                sample.name, reformatted_forward, reformatted_reverse
            )
            # Run repair.sh to repair the reads
            out, err, _ = bbtools.repair_reads(
                forward_in=reformatted_forward,
                reverse_in=reformatted_reverse,
                forward_out=repair_forward,
                reverse_out=repair_reverse,
                returncmd=True
            )
            logger.debug('repair.sh completed for sample: %s', sample.name)

            # Write the output and error messages to the log file
            logger.debug(
                'Writing repair.sh output to log file for sample: %s',
                sample.name
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

            # Update the fastq_files attribute to point to the repaired files
            logger.debug(
                'Updating fastq_files attribute for sample: %s', sample.name
            )
            sample.general.fastq_files = (
                [repair_forward, repair_reverse]
                if repair_reverse else [repair_forward]
            )
        except CalledProcessError as exc:
            # Log an error if the repair process fails
            logger.error(
                'repair.sh failed for sample: %s with error: %s',
                sample.name, str(exc)
            )
            logger.error(traceback.format_exc())
            log_fastq_error(
                log_file=log_file,
                logger=logger,
                sample=sample
            )
            error(
                logger=logger,
                message='fastq_error',
                sample=sample
            )
    else:
        # Log an error if the reformatted files do not exist
        logger.error(
            'Reformatted files do not exist for sample: %s', sample.name
        )
        log_fastq_error(
            log_file=log_file,
            logger=logger,
            sample=sample
        )
        error(
            logger=logger,
            message='fastq_error',
            sample=sample
        )


def log_fastq_error(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Log an error message for FASTQ file issues.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger
        sample (CustomBox): Metadata sample object.
    """
    logger.debug('Logging FASTQ error for sample: %s', sample.name)

    # Create an error message indicating issues with the FASTQ files
    message = (
        f'An error was detected in the FASTQ files for sample {sample.name}. '
        'These files will not be processed further'
    )
    logger.error(message)

    # Write the error message to the log file and sample-specific logs
    logger.debug(
        'Writing FASTQ error message to log file for sample: %s', sample.name
    )
    write_to_log_file(
        out=message,
        err=message,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )


def error(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    message: str,
    sample: CustomBox
) -> None:
    """
    Check to see if the run CustomBox exists. If so, update the
    run.status to reflect the error.

    Args:
        logger (logging.Logger): Logger object.
        message (str): Error message to add to the sample.run.status attribute.
        sample (CustomBox): Metadata sample object.
    """
    logger.debug('Setting error status for sample: %s', sample.name)

    # Set the .fastq_files attribute to an empty list
    sample.general.fastq_files = []
    logger.debug('Cleared fastq_files for sample: %s', sample.name)

    # Ensure that the run attribute exists and is a CustomBox
    if not hasattr(sample, 'run'):
        sample.run = CustomBox()
        logger.debug('Created run attribute for sample: %s', sample.name)

    # Set the status attribute in the run CustomBox
    sample.run.status = message
    logger.debug(
        'Set run.status to "%s" for sample: %s', message, sample.name
    )
