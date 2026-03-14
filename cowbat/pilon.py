#!/usr/env/python3

"""
Run pilon misassembly fixing on the assemblies
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


def pilon(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> None:
    """
    Run pilon to fix any misassemblies in the contigs - will look for SNPs
    and indels.

    Arguments:
    - `log_file` should be the path to the log file.
    - `logger` should be the logger object.
    - `metadata` should contain the sample data with necessary
        attributes.
    - `threads` should be set to the number of threads to use.
    """
    logger.info('Improving quality of assembly with pilon')

    # Process each sample
    for sample in metadata:
        if _is_valid_pilon_sample(logger=logger, sample=sample):
            logger.debug(
                "Initializing pilon for sample: %s", sample.name)
            sample = _initialize_sample_pilon(
                logger=logger,
                sample=sample,
                threads=threads
            )
            logger.debug(
                "Submitting pilon task for sample: %s", sample.name
            )
            _run_pilon(
                log_file=log_file,
                logger=logger,
                sample=sample
            )

    return metadata


def _is_valid_pilon_sample(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> bool:
    """
    Check if the sample is valid for running pilon.

    Args:
        logger: The logger object.
        sample: The sample object to check.

    Returns:
        bool: True if the sample is valid, False otherwise.
    """
    logger.debug("Checking if sample is valid for pilon: %s", sample.name)

    # Check if the best assembly file is not 'NA'
    if sample.general.best_assembly_file == 'NA':
        logger.debug(
            "Sample %s is not valid for pilon: best_assembly_file is 'NA'",
            sample.name)
        return False

    # Check if the necessary attributes are present
    if not hasattr(
            sample,
            'quast') or not hasattr(
            sample.quast,
            'sorted_bam'):
        logger.debug(
            "Sample %s is not valid for pilon: missing quast outputs or "
            "sorted_bam file", sample.name
        )
        return False
    if not hasattr(
            sample,
            'general') or not hasattr(
            sample.general,
            'assembly_file'):
        logger.debug(
            "Sample %s is not valid for pilon: missing general or "
            "assembly_file", sample.name
        )
        return False
    return True


def _initialize_sample_pilon(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    sample: CustomBox,
    threads: int
) -> None:
    """
    Initialize the pilon attributes for a sample.

    Args:
        logger: The logger object.
        sample: The sample object to initialize attributes for.
        threads: The number of threads to use in the analyses
    """
    logger.debug(
        "Initializing pilon attributes for sample: %s",
        sample.name)
    # Initialize the pilon attribute
    sample.pilon = CustomBox()

    # Set the contigs file to the assembly file
    sample.general.contigs_file = sample.general.assembly_file

    # Set the output directory for pilon
    sample.pilon.out_dir = os.path.join(sample.quast.output_dir, 'pilon')
    os.makedirs(sample.pilon.out_dir, exist_ok=True)

    # Construct the pilon command using f-strings
    sample.pilon.cmd = (
        f'pilon --genome {sample.general.contigs_file} '
        f'--bam {sample.quast.sorted_bam} --fix bases '
        f'--threads {threads} --out_dir {sample.pilon.out_dir} '
        f'--changes --mindepth 0.25'
    )
    logger.debug(
        "Pilon command for sample %s: %s",
        sample.name,
        sample.pilon.cmd)

    return sample


def _run_pilon(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Run pilon for a sample.

    Args:
        log_file: The name and path of the log file
        logger: The logger object.
        sample: The sample object to run pilon on.
    """
    logger.debug("Running pilon for sample: %s", sample.name)
    # Check if the contigs file already exists
    if not os.path.isfile(sample.general.contigs_file):
        logger.debug(
            "Contigs file does not exist for sample: %s",
            sample.name)
        _execute_pilon_command(
            log_file=log_file,
            logger=logger,
            sample=sample
        )


def _execute_pilon_command(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Execute the pilon command for a sample.

    Args:
        log_file: The name and path of the log file
        logger: The logger object.
        sample: The sample object to run the command on.
    """
    try:
        logger.debug(
            "Executing pilon command for sample: %s",
            sample.name)
        # Run the pilon command
        out, err = run_subprocess(sample.pilon.cmd)
        # Log the output and errors
        _log_pilon_output(
            err=err,
            log_file=log_file,
            logger=logger,
            out=out,
            sample=sample
        )
    except RuntimeError as exc:
        logger.error(
            "Failed to run pilon command for sample %s: %s",
            sample.name, exc
        )


def _log_pilon_output(
    *,  # Enforce keyword arguments
    err: str,
    log_file: str,
    logger: logging.Logger,
    out: str,
    sample: CustomBox
) -> None:
    """
    Log the output and errors from the pilon command.

    Args:
        log_file: The name and path of the log file
        logger: The logger object.
        sample: The sample object to log information for.
        out: The standard output from the pilon command.
        err: The standard error from the pilon command.
    """
    logger.debug("Logging pilon output for sample: %s", sample.name)
    write_to_log_file(
        out=f'{sample.pilon.cmd}\n{out}',
        err=err,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err,
    )
