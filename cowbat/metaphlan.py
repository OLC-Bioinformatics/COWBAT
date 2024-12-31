#!/usr/bin/env python3

"""
Run metaphlan analyses on the samples
"""

# Standard imports
import logging
import os
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_metadata_to_file,
)
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.methods import write_to_log_files

__author__ = 'adamkoziol'


def run_metaphlan_analyses(
    *,  # Enforce keyword arguments
    error_logger: logging.Logger,
    file_format: str,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    report_path: str,
    sequence_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Run metaphlan analyses on the samples

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        file_format (str): Format of the sequence files.
            Either 'fasta' or 'fastq'.
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
        RuntimeError: If there is an issue with the metaphlan analyses.
    """
    # Create the report path
    os.makedirs(report_path, exist_ok=True)

    # Iterate through the samples
    for sample in metadata:
        # Skip samples without sequence data when the file format is 'fasta'
        if (
            sample.general.best_assembly_file == 'NA' and
            file_format == 'fasta'
        ):
            continue

        # Set the metaphlan attributes
        sample = _set_metadata_attributes(
            file_format=file_format,
            sample=sample,
            sequence_path=sequence_path
        )

        # Set the metaphlan command
        sample = _set_metaphlan_command(
            file_format=file_format,
            sample=sample,
            threads=threads
        )

        logger.debug(
            'Metaphlan command: %s for sample %s',
            sample.commands.metaphlan[file_format], sample.name
        )

        # Run the metaphlan command if the report does not exist
        if os.path.isfile(sample.metaphlan[file_format].report):
            continue

        out, err = run_subprocess(
            command=sample.commands.metaphlan[file_format]
        )

        # Write the command and the outputs to the log files
        write_to_log_files(
            command=sample.commands.metaphlan[file_format],
            err=err,
            log_file=log_file,
            logger=logger,
            out=out,
            program='metaphlan',
            sample=sample
        )

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        logger=logger,
        metadata=metadata
    )

    # Return the updated metadata
    return metadata


def _set_metadata_attributes(
    *,  # Enforce keyword arguments
    file_format: str,
    sequence_path: str,
    sample: CustomBox
) -> CustomBox:
    """
    Set the metaphlan attributes for the sample

    Args:
        file_format (str): The format of the sequence files.
        sequence_path (str): The path to the sequence files.
        sample (CustomBox): The sample object.

    Returns:
        CustomBox: The updated sample object.
    """
    # Create the metaphlan attribute
    if not sample.key_exists('metaphlan'):
        sample.metaphlan = CustomBox()

    # Create the metaphlan file_format object
    if not sample.metaphlan.key_exists(file_format):
        sample.metaphlan[file_format] = CustomBox()

    # Create the metaphlan output directory
    sample.metaphlan[file_format].report_dir = os.path.join(
        sample.general.output_directory, f'metaphlan_{file_format}'
    )
    os.makedirs(sample.metaphlan[file_format].report_dir, exist_ok=True)

    # Set the name of the metaphlan report
    sample.metaphlan[file_format].report = os.path.join(
        sample.metaphlan[file_format].report_dir,
        f'{sample.name}.{file_format}.metaphlan.txt'
    )

    # Set the name of the bowtie2 output file for the FASTQ files
    if file_format == 'fastq':
        sample.metaphlan[file_format].bowtie2out = os.path.join(
            sample.metaphlan[file_format].report_dir,
            f'{sample.name}.{file_format}.bowtie2.bz2'
        )
    elif file_format == 'fasta':
        # Set the path of the BestAssemblies folder
        assembly_path = os.path.join(sequence_path, 'BestAssemblies')

        # Set the name of the bowtie2 output file for the assembly
        sample.metaphlan[file_format].bowtie2out = os.path.join(
            assembly_path,
            f'{sample.name}.{file_format}.bowtie2out.txt'
        )

    # Create the file format specific attributes for the commands
    sample.commands.metaphlan = CustomBox()

    return sample


def _set_metaphlan_command(
    *,  # Enforce keyword arguments
    file_format: str,
    sample: CustomBox,
    threads: int
) -> CustomBox:
    """
    Set the metaphlan command for the sample

    Args:
        file_format (str): The format of the sequence files.
        sample (CustomBox): The sample object.
        threads (int): The number of threads to use for processing.

    Returns:
        CustomBox: The updated sample object.
    """
    # If the bowtie2 output file exists, run the metaphlan command using the
    # bowtie2 output file
    if os.path.isfile(sample.metaphlan[file_format].bowtie2out):
        sample.commands.metaphlan[file_format] = (
            f'metaphlan {sample.metaphlan[file_format].bowtie2out} '
            f'--nproc {threads} --input_type bowtie2out '
            f'-o {sample.metaphlan[file_format].report}'
        )
        return sample

    # If the bowtie2 output file does not exist, run the metaphlan command for
    # the appropriate file format
    if file_format == 'fastq':
        sample.commands.metaphlan[file_format] = (
            f'metaphlan '
            f'{",".join(sample.general.trimmed_corrected_fastq_files)} '
            f'--input_type fastq '
            f'--bowtie2out {sample.metaphlan[file_format].bowtie2out} '
            f'--nproc {threads} '
            f'-o {sample.metaphlan[file_format].report}'
        )
    elif file_format == 'fasta':
        sample.commands.metaphlan[file_format] = (
            f'metaphlan {sample.general.best_assembly_file} '
            f'--input_type fasta '
            f'--nproc {threads} '
            f'-o {sample.metaphlan[file_format].report}'
        )
    else:
        raise ValueError(f'Invalid format: {file_format}')

    return sample
