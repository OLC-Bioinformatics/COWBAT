#!/usr/bin/env python3

"""
Assembly evaluation using various tools including Bowtie2, Samtools, Quast,
and Qualimap.
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


def prepare_sample_for_build(
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Prepare the sample for the bowtie2-build command.

    Args:
        logger (logging.Logger): Logger for recording information.
        sample (CustomBox): Metadata sample object.
    """
    sample.quast = CustomBox()
    sample.quast.base_name = os.path.splitext(sample.general.assembly_file)[0]
    sample.quast.build_command = (
        f'bowtie2-build {sample.general.assembly_file} '
        f'{sample.quast.base_name}'
    )
    logger.debug("Constructed build command: %s", sample.quast.build_command)


def index_files_exist(sample: CustomBox) -> bool:
    """
    Check if the index files already exist.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        bool: True if index files exist, False otherwise.
    """
    return os.path.isfile(f'{sample.quast.base_name}.1.bt2')


def run_build_command(
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Run the bowtie2-build command for the sample.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        sample (CustomBox): Metadata sample object.

    Raises:
        RuntimeError: If the subprocess command fails.
    """
    try:
        out, err = run_subprocess(sample.quast.build_command)
        logger.debug("Command output: %s", out)
        logger.debug("Command error: %s", err)

        # Write the appropriate information to the log_file
        write_to_log_file(
            out=f'{sample.quast.build_command}\n{out}',
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )
    except RuntimeError as exc:
        logger.error(
            "Failed to run build command for sample %s: %s",
            sample.name, exc
        )


def bowtie_build(
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Use bowtie2-build to index each target file.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata sample objects.

    Raises:
        FileNotFoundError: If any required file paths are not found.
        RuntimeError: If the subprocess command fails.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Preparing targets for reference mapping')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)
        prepare_sample_for_build(
            logger=logger,
            sample=sample
        )

        if not index_files_exist(sample):
            logger.info(
                "Index files not found for sample %s. Running bowtie2-build",
                sample.name
            )
            run_build_command(
                log_file=log_file,
                logger=logger,
                sample=sample
            )
        else:
            logger.info(
                "Index files already exist for sample %s, skipping build",
                sample.name
            )

    return metadata


def prepare_output_directory(sample: CustomBox) -> None:
    """
    Prepare the output directory for the sample.

    Args:
        sample (CustomBox): Metadata sample object.
    """
    sample.quast.output_dir = os.path.join(
        sample.general.output_directory, 'quast'
    )
    os.makedirs(sample.quast.output_dir, exist_ok=True)
    sample.quast.sorted_bam = os.path.join(
        sample.quast.output_dir, f'{sample.name}_sorted.bam'
    )


def construct_mapping_command(
    logger: logging.Logger,
    sample: CustomBox,
    threads: int
) -> None:
    """
    Construct the Bowtie2 mapping command for the sample.

    Args:
        logger (logging.Logger): Logger for recording information.
        sample (CustomBox): Metadata sample object.
        threads (int): Number of threads to use.
    """
    sample.quast.map_command = f'bowtie2 -x {sample.quast.base_name}'

    # Add FASTQ files to the command
    if len(sample.general.trimmed_corrected_fastq_files) == 1:
        sample.quast.map_command += (
            f' -U {sample.general.trimmed_corrected_fastq_files[0]}'
        )
    else:
        sample.quast.map_command += (
            f' -1 {sample.general.trimmed_corrected_fastq_files[0]} '
            f'-2 {sample.general.trimmed_corrected_fastq_files[1]}'
        )

    # Add Samtools commands to convert and sort the BAM file
    sample.quast.map_command += (
        f' -p {threads} -X 1000 | '
        f'samtools view -@ {threads} -h -F 4 -bT '
        f'{sample.general.assembly_file} - | '
        f'samtools sort - -@ {threads} '
        f'-o {sample.quast.sorted_bam}'
    )
    logger.debug("Constructed map command: %s", sample.quast.map_command)


def sorted_bam_exists(sample: CustomBox) -> bool:
    """
    Check if the sorted BAM file already exists.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        bool: True if the sorted BAM file exists, False otherwise.
    """
    return os.path.isfile(sample.quast.sorted_bam)


def run_mapping_command(
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Run the Bowtie2 mapping command for the sample.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        sample (CustomBox): Metadata sample object.

    Raises:
        RuntimeError: If the subprocess command fails.
    """
    try:
        out, err = run_subprocess(sample.quast.map_command)
        logger.debug("Command output: %s", out)
        logger.debug("Command error: %s", err)

        # Write the appropriate information to the log_file
        write_to_log_file(
            out=f'{sample.quast.map_command}\n{out}',
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except RuntimeError as exc:
        logger.error(
            "Failed to run mapping command for sample %s: %s",
            sample.name, exc
        )


def bowtie_run(
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Map the FASTQ reads against the appropriate target file using Bowtie2
    and Samtools.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata sample objects.
        threads (int): Number of threads to use.

    Raises:
        FileNotFoundError: If any required file paths are not found.
        RuntimeError: If the subprocess command fails.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Performing reference mapping for quality evaluation')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        if sample.general.best_assembly_file == "NA":
            logger.warning(
                "Skipping sample %s as best assembly file is 'NA'",
                sample.name
            )
            continue

        prepare_output_directory(
            sample=sample
        )
        construct_mapping_command(
            logger=logger,
            sample=sample,
            threads=threads
        )

        if not sorted_bam_exists(sample):
            logger.info(
                "Sorted BAM file not found for sample %s, running command",
                sample.name
            )
            run_mapping_command(
                log_file=log_file,
                logger=logger,
                sample=sample
            )
        else:
            logger.info(
                "Sorted BAM file already exists for sample %s, skipping",
                sample.name
            )

    return metadata
