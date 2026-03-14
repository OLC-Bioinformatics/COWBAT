#!/usr/env/python3

"""
Annotate genomes with Bakta
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
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def bakta(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    sequence_path: str,
) -> List[CustomBox]:
    """
    Annotation of genomes with Bakta

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        sequence_path (str): Path to the sequence files.

    Returns:
        List[CustomBox]: Updated metadata after all processing steps.

    Raises:
        RuntimeError: If there is an issue with the Bakta analyses.
    """

    # Create the necessary attributes
    metadata = _create_attributes(
        metadata=metadata,
        reference_file_path=reference_file_path,
        sequence_path=sequence_path
    )

    for sample in metadata:
        # Run Bakta analyses
        if sample.general.best_assembly_file == 'NA':
            continue

        logger.debug('Running Bakta command %s', sample.commands.bakta)

        # Run Bakta
        out, err = run_subprocess(
            command=sample.commands.bakta
        )

        # Write the command to the log file
        write_to_log_file(
            out=sample.commands.bakta,
            err=sample.commands.bakta,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )

        # Write the stdout and stderr to the log file
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )

    return metadata


def _create_attributes(
    *,  # Enforce keyword arguments
    metadata: List[CustomBox],
    reference_file_path: str,
    sequence_path: str
) -> List[CustomBox]:
    """
    Create the necessary attributes for Bakta analyses

    Args:
        metadata (List[CustomBox]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        sequence_path (str): Path to the sequence files.

    Returns:
        List[CustomBox]: Updated metadata
    """
    # Set the bakta-specific reference file path
    bakta_reference_file_path = os.path.join(
        reference_file_path,
        'bakta',
        'db'
    )

    # Create the necessary attributes
    for sample in metadata:
        sample.bakta = CustomBox()
        sample.bakta.output_dir = os.path.join(
            sample.general.output_directory,
            'bakta'
        )

        # Set the docker-specific paths
        sample.bakta.docker_assembly_file = os.path.join(
            '/data',
            'BestAssemblies',
            sample.name + '.fasta'
        )
        sample.bakta.docker_output_dir = os.path.join(
            '/data',
            sample.name,
            'bakta'
        )

        # Create the output directory
        os.makedirs(sample.bakta.output_dir, exist_ok=True)

        # Create the system call
        sample.commands.bakta = (
            f'docker run -v {bakta_reference_file_path}:/db '
            f'-v {sequence_path}:/data '
            '--rm --entrypoint /bin/bash oschwengers/bakta:latest '
            f'-c "bakta '
            f'-o {sample.bakta.docker_output_dir} '
            f'{sample.bakta.docker_assembly_file} '
            f'--force --complete --db /db"'
        )
        print(sample.commands.bakta)

    return metadata
