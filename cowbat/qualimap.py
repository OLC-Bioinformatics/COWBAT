#!/usr/env/python3

"""
Qualimap functions
"""

# Standard imports
from concurrent.futures import as_completed, ThreadPoolExecutor
import logging
import os
from typing import Dict, List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox


def qualimapper(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Create threads and commands for performing reference mapping for
    qualimap analyses.

    This method initializes threads and enqueues tasks for qualimap
    analysis. It waits for all tasks to be completed.

    Arguments:
    - `log_file` should be the path to the log file.
    - `logger` should be the logger object.
    - `metadata` should contain the sample data with necessary
        attributes.
    - `threads` should be set to the number of threads to use.

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of sample metadata.
        threads (int): Number of threads to use for processing.

    Returns:
        List[CustomBox]: Updated metadata after qualimap analyses.
    """
    logger.info('Running qualimap on samples')

    futures = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for sample in metadata:
            sample = _initialize_sample_qualimap(sample=sample)
            if sample.general.best_assembly_file != "NA":
                future = executor.submit(_run_qualimap, log_file, sample)
                futures.append(future)

    updated_metadata = []
    for future in as_completed(futures):
        try:
            result = future.result()
            updated_metadata.append(result)
        except (IOError, OSError) as file_error:
            logger.error("File operation error: %s", file_error)
        except RuntimeError as runtime_error:
            logger.error(
                "Runtime error during qualimap analysis: %s", runtime_error
            )

    return updated_metadata


def _initialize_sample_qualimap(*, sample: CustomBox) -> None:
    """
    Initialize the qualimap attributes for a sample.

    Args:
        sample: The sample object to initialize attributes for.
    """
    # Create and populate the qualimap attribute
    sample.qualimap = CustomBox()
    sample.qualimap.output_dir = os.path.join(
        sample.general.output_directory, 'qualimap'
    )
    os.makedirs(sample.qualimap.output_dir, exist_ok=True)
    sample.qualimap.report_file = os.path.join(
        sample.qualimap.output_dir, 'genome_results.txt'
    )

    # Initialize dictionaries to store qualimap results
    sample.qualimap.length = {}
    sample.qualimap.bases = {}
    sample.qualimap.coverage = {}
    sample.qualimap.std_dev = {}

    return sample


def _run_qualimap(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Run qualimap for a sample.

    Args:
        log_file: The name and path of the log file
        logger: The logger object to record information.
        sample: The sample object to run qualimap on.
    """
    if not _is_valid_sample(sample=sample):
        logger.warning(
            "Invalid sample %s, skipping qualimap run", sample.name
        )
        return

    qualimap_call = _construct_qualimap_command(sample=sample)
    sample.commands.qualimap = qualimap_call

    logger.debug(
        'Qualimap call for sample %s: %s',
        sample.name, sample.commands.qualimap
    )

    if not os.path.isfile(sample.qualimap.report_file):
        _execute_qualimap_command(
            log_file=log_file,
            logger=logger,
            sample=sample
        )

    return sample


def _is_valid_sample(*, sample: CustomBox) -> bool:
    """
    Check if the sample is valid for running qualimap.

    Args:
        sample: The sample object to check.

    Returns:
        bool: True if the sample is valid, False otherwise.
    """
    if sample.general.best_assembly_file == "NA":
        return False
    if not hasattr(
            sample,
            'quast') or not hasattr(
            sample.quast,
            'sorted_bam'):
        return False
    if not hasattr(
            sample,
            'qualimap') or not hasattr(
            sample.qualimap,
            'output_dir'):
        return False
    return True


def _construct_qualimap_command(*, sample: CustomBox) -> str:
    """
    Construct the qualimap command for a sample.

    Args:
        sample: The sample object to construct the command for.

    Returns:
        str: The constructed qualimap command.
    """
    return (
        f'qualimap bamqc -bam {sample.quast.sorted_bam} '
        f'-outdir {sample.qualimap.output_dir}'
    )


def _execute_qualimap_command(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Execute the qualimap command for a sample.

    Args:
        log_file: The name and path of the log file
        logger: The logger object to record information.
        sample: The sample object to run the command on.
    """
    try:
        out, err = run_subprocess(sample.commands.qualimap)
        _log_qualimap_output(
            err=err,
            log_file=log_file,
            out=out,
            sample=sample
        )
    except RuntimeError as exc:
        logger.error(
            "Failed to run qualimap command for sample %s: %s",
            sample.name, exc
        )


def _log_qualimap_output(
    *,  # Enforce keyword arguments
    err: str,
    log_file: str,
    out: str,
    sample: CustomBox
) -> None:
    """
    Log the output and errors from the qualimap command.

    Args:
        log_file: The name and path of the log file
        sample: The sample object to log information for.
        out: The standard output from the qualimap command.
        err: The standard error from the qualimap command.
    """
    write_to_log_file(
        out=f'{sample.commands.qualimap}\n{out}',
        err=err,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err,
    )


def parse_qualimap_report(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> None:
    """
    Parse the qualimap report.

    This method iterates over the metadata samples, reads the qualimap
    report, and populates the metadata object with the extracted key-value
    pairs.

    Arguments:
    - `logger` should be the logger
    - `metadata` should contain the sample data with necessary
        attributes.
    """
    logger.info('Parsing Qualimap reports')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        # Check if the sample is valid for parsing
        if not _is_valid_qualimap_sample(sample=sample):
            logger.warning(
                "Skipping sample %s as best assembly file is 'NA' or "
                "report file not found", sample.name
            )
            continue

        # Parse the report file and extract key-value pairs
        qualimap_dict = _parse_report_file(
            logger=logger,
            sample=sample
        )
        # Update the sample object with the extracted key-value pairs
        sample = _update_sample_with_qualimap_dict(
            qualimap_dict=qualimap_dict,
            sample=sample
        )

    return metadata


def _is_valid_qualimap_sample(*, sample: CustomBox) -> bool:
    """
    Check if the sample is valid for parsing qualimap report.

    Args:
        sample: The sample object to check.

    Returns:
        bool: True if the sample is valid, False otherwise.
    """
    # Check if the best assembly file is not 'NA'
    if sample.general.best_assembly_file == "NA":
        return False
    # Check if the qualimap report file exists
    if not os.path.isfile(sample.qualimap.report_file):
        return False
    return True


def _parse_report_file(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> dict:
    """
    Parse the qualimap report file and extract key-value pairs.

    Args:
        logger: The logger object to record information.
        sample: The sample object to parse the report for.

    Returns:
        A dictionary containing the extracted key-value pairs.
    """
    qualimap_dict = {}
    try:
        with open(sample.qualimap.report_file, encoding='utf-8') as report:
            for line in report:
                # Sanitize the keys and values using qualimap_analyze
                key, value = qualimap_analyze(line=line)

                # If the keys and values exist, enter them into the dict
                if (key, value) != (None, None):
                    # Only keep two decimal places for float values
                    if isinstance(value, float):
                        value = float(f'{value:.2f}')
                    qualimap_dict[key] = value

                # Parse the coverage per contig section
                if 'Coverage per contig' in line:
                    _parse_contig_coverage(
                        logger=logger,
                        report=report,
                        sample=sample
                    )
    except (IOError, FileNotFoundError) as exc:
        logger.error(
            "Error reading qualimap report file for sample %s: %s",
            sample.name, exc
        )
    return qualimap_dict


def _parse_contig_coverage(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    report: str,
    sample: CustomBox
) -> None:
    """
    Parse the coverage per contig section of the qualimap report.

    Args:
        logger: The logger object to record information.
        report: The file object for the qualimap report.
        sample: The sample object to populate with extracted data.
    """
    for contig_line in report:
        try:
            # Extract the contig coverage details
            _, name, length, bases, coverage, std_dev = (
                contig_line.rstrip().split('\t')
            )

            # Update the sample object with the extracted details
            sample.qualimap.length[name] = length
            sample.qualimap.bases[name] = bases
            sample.qualimap.coverage[name] = coverage
            sample.qualimap.std_dev[name] = std_dev
        except ValueError:
            logger.debug(
                "Skipping malformed line in contig coverage section: %s",
                contig_line
            )


def _update_sample_with_qualimap_dict(
    *,  # Enforce keyword arguments
    qualimap_dict: Dict,
    sample: CustomBox,
) -> None:
    """
    Update the sample object with the extracted key-value pairs.

    Args:
        sample: The sample object to update.
        qualimap_dict: The dictionary containing the extracted key-value
            pairs.
    """
    if qualimap_dict:
        for attribute, value in qualimap_dict.items():
            # Remove the 'X' from the depth values e.g. 40.238X
            setattr(sample.qualimap, attribute, value.rstrip('X'))

    return sample


def qualimap_analyze(*, line: str):
    """
    Analyze a line from the qualimap report.

    Args:
        line: A line from the qualimap report.

    Returns:
        A tuple containing the sanitized key and value.
    """
    # Check if the line contains ' = '
    if ' = ' in line:
        key, value = line.split(' = ', 1)

        # Sanitize the key
        key = (
            key.replace('number of ', "")
            .replace("'", "")
            .title()
            .replace(" ", "")
        )

        # Sanitize the value
        value = value.replace(",", "").replace(" ", "").rstrip()
    else:
        # Set the keys and values to None if ' = ' is not found
        key, value = None, None

    return key, value
