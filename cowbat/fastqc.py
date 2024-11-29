#!/usr/bin/env python3

"""
FastQC functions
"""

# Standard imports
from concurrent.futures import as_completed, ThreadPoolExecutor
from glob import glob
import logging
import os
import shutil
from queue import Queue
from typing import List, Tuple

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file,
)
from olctools.accessoryFunctions.metadata import CustomBox

qc_queue = Queue()
__author__ = 'adamkoziol'


def fastqc_threader(
    *,  # Enforce keyword arguments
    level: str,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Run quality control on FASTQ files using FastQC.

    Args:
        level (str): The level of processing (e.g., 'Trimmed', 'merged').
        log_file (str): Path to the log file.
        metadata (List[CustomBox]): List of metadata sample objects.
        threads (int): Number of threads to use.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Running quality control on %s fastq files', level)

    # Create and start threads for each FASTQ file in the list
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                fastqc,
                sample=sample,
                level=level,
                log_file=log_file,
                threads=threads) for sample in metadata if isinstance(
                sample.general.fastq_files,
                list)]
        logger.info('Submitted %d tasks to the executor', len(futures))

    # Collect results from futures
    updated_metadata = []
    for future in as_completed(futures):
        result = future.result()
        updated_metadata.append(result)
        logger.info('Task completed for sample: %s', result.name)

    # Wait on the queue until everything has been processed
    logger.info('Waiting for the queue to be processed')
    while not qc_queue.empty():
        qc_queue.get()
        qc_queue.task_done()
    logger.info('Queue processing completed')

    return updated_metadata


def fastqc(
    *,  # Enforce keyword arguments
    level: str,
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox,
    threads: int
) -> CustomBox:
    """
    Run FastQC on the given sample.

    Args:
        level (str): The level of processing (e.g., 'Trimmed', 'merged').
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording errors.
        sample (CustomBox): Metadata sample object.
        threads (int): Number of threads to use.

    Returns:
        CustomBox: Updated metadata sample object.
    """
    logger.debug('Preparing FastQC call for sample: %s', sample.name)

    # Prepare the FastQC system call and reads call
    fastqc_call, fastqc_reads = prepare_fastqc_call(
        level=level,
        logger=logger,
        sample=sample,
        threads=threads
    )

    # If a FastQC call is prepared, add it to the queue
    if fastqc_call:
        sample.commands.fastqc = fastqc_call
        setattr(sample.commands, f'fastqc_{level.lower()}', fastqc_call)
        qc_queue.put(
            (sample, fastqc_call, fastqc_reads, log_file, level, logger)
        )
        logger.debug('Added FastQC call to queue for sample: %s', sample.name)

    # Process the FastQC queue
    while not qc_queue.empty():
        process_fastqc_queue()

    return sample


def prepare_fastqc_call(
    *,  # Enforce keyword arguments
    level: str,
    logger: logging.Logger,
    sample: CustomBox,
    threads: int
) -> Tuple[str, str]:
    """
    Prepare the FastQC system call based on the processing level.

    Args:
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        logger (logging.Logger): Logger for recording errors.
        sample (CustomBox): Metadata sample object.
        threads (int): Number of threads to use.

    Returns:
        tuple: FastQC system call and reads call.
    """
    logger.debug(
        'Preparing FastQC call for level: %s, sample: %s',
        level,
        sample.name)

    reader = 'cat'
    fastq_files = None

    # Determine the reader and FASTQ files based on the processing level
    if level == 'trimmed':
        reader, fastq_files = get_reader_and_files(
            attr='trimmed_fastq_files',
            logger=logger,
            sample=sample
        )
    elif level == 'trimmed_corrected':
        reader, fastq_files = get_reader_and_files(
            attr='trimmed_corrected_fastq_files',
            logger=logger,
            sample=sample
        )
    elif level == 'normalised':
        reader, fastq_files = get_reader_and_files(
            attr='normalised_reads',
            logger=logger,
            sample=sample
        )
    elif level == 'merged':
        reader, fastq_files = get_reader_and_files(
            attr='merged_reads',
            logger=logger,
            sample=sample,
            single_file=True
        )
    else:
        reader, fastq_files = get_reader_and_files(
            attr='fastq_files',
            logger=logger,
            sample=sample
        )

    # If no valid FASTQ files are found, return empty strings
    if not isinstance(fastq_files, list):
        logger.debug(
            'No valid %s FASTQ files found for sample: %s',
            level, sample.name
        )
        return '', ''

    # Set the output directory for FastQC results
    out_dir = os.path.join(sample.general.output_directory, 'fastqc', level)
    os.makedirs(out_dir, exist_ok=True)

    # Prepare the FastQC system call and reads call based on the number of
    # FASTQ files
    if len(fastq_files) == 2:
        fastqc_call = (
            f'{reader} {fastq_files[0]} {fastq_files[1]} | fastqc -q -t '
            f'{threads} stdin -o {out_dir}'
        )
        fastqc_reads = (
            f"fastqc {fastq_files[0]} {fastq_files[1]} -q -o {out_dir} -t "
            f"{threads}"
        )
    elif len(fastq_files) == 1:
        fastqc_call = (
            f'{reader} {fastq_files[0]} | fastqc -q -t {threads} stdin -o '
            f'{out_dir}'
        )
        fastqc_reads = (
            f"fastqc {fastq_files[0]} -q -o {out_dir} -t {threads}"
        )
    else:
        fastqc_call = ''
        fastqc_reads = ''

    logger.debug(
        'Prepared FastQC call for %s reads for sample: %s',
        level, sample.name
    )
    return fastqc_call, fastqc_reads


def get_reader_and_files(
    *,  # Enforce keyword arguments
    attr: str,
    logger: logging.Logger,
    sample: CustomBox,
    single_file: bool = False
) -> Tuple[str, List[str]]:
    """
    Get the appropriate reader and FASTQ files based on the attribute.

    Args:
        attr (str): Attribute name to get the FASTQ files.
        logger (logging.Logger): Logger for recording errors.
        sample (CustomBox): Metadata sample object.
        single_file (bool): Whether to expect a single file or a list of files.

    Returns:
        tuple: Reader command and FASTQ files.
    """
    try:
        # Get the FASTQ files from the specified attribute
        fastq_files = getattr(sample.general, attr)

        if single_file:
            fastq_files = [fastq_files]

        # Determine the appropriate reader based on the file extension
        if '.gz' in fastq_files[0]:
            reader = 'gunzip --to-stdout'
        elif '.bz2' in fastq_files[0]:
            reader = 'bunzip2 --stdout'
        else:
            reader = 'cat'
    except AttributeError:
        logger.debug(
            'AttributeError: %s not found in sample: %s',
            attr,
            sample.name
        )
        reader = 'cat'
        fastq_files = []

    logger.debug(
        'Reader: %s, FASTQ files: %s for sample: %s',
        reader,
        fastq_files,
        sample.name
    )
    return reader, fastq_files


def process_fastqc_queue() -> None:
    """
    Process the FastQC queue.
    """
    while not qc_queue.empty():
        # Get the next item from the queue
        sample, system_call, fastqc_reads, log_file, level, logger = \
            qc_queue.get()
        output_dir = os.path.join(sample.general.output_directory, 'fastqc')

        try:
            # Check if the FastQC output HTML file already exists
            _ = glob(os.path.join(output_dir, '*.html'))[0]
        except IndexError:
            # If the output file does not exist, create the output directory
            os.makedirs(output_dir, exist_ok=True)

            # Run the FastQC system calls and log the output
            run_fastqc(
                fastqc_reads=fastqc_reads,
                level=level,
                log_file=log_file,
                logger=logger,
                output_dir=output_dir,
                sample=sample,
                system_call=system_call,
            )
            # Rename the FastQC output files
            rename_fastqc_outputs(
                level=level,
                logger=logger,
                output_dir=output_dir,
                sample=sample
            )

        # Mark the task as done
        qc_queue.task_done()


def run_fastqc(
    *,  # Enforce keyword arguments
    fastqc_reads: str,
    log_file: str,
    logger: logging.Logger,
    level: str,
    output_dir: str,
    sample: CustomBox,
    system_call: str,
) -> Tuple[str, str]:
    """
    Run the FastQC system calls and log the output.

    Args:
        system_call (str): FastQC system call.
        fastqc_reads (str): FastQC reads call.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording errors.
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        output_dir (str): Output directory.
        sample (CustomBox): Metadata sample object.

    Returns:
        tuple: Standard output and error strings.
    """
    logger.debug(
        'Running FastQC system call for %s reads for sample: %s',
        level, sample.name
    )

    # Set the name of the output file
    fastqc_html = os.path.join(output_dir, f'{sample.name}_fastqc.html')

    # Only run the system call if the FastQC outputs don't already exist
    if os.path.isfile(fastqc_html):
        logger.debug(
            'FastQC output already exists for sample: %s, skipping FastQC',
            sample.name
        )
        return "", ""

    out_str, err_str = '', ''

    # Run the first FastQC system call
    out, err = _run_command(
        command=system_call,
        logger=logger
    )
    out_str += out
    err_str += err

    logger.debug(
        'Running FastQC reads call for %s reads for sample: %s',
        level, sample.name
    )

    # Run the second FastQC reads call
    out, err = _run_command(
        command=fastqc_reads,
        logger=logger
    )
    out_str += ' ' + out
    err_str += ' ' + err

    # Write the system call and reads call to the log file
    _log_fastqc_output(
        err_str=err_str,
        fastqc_reads=fastqc_reads,
        log_file=log_file,
        logger=logger,
        out_str=out_str,
        sample=sample,
        system_call=system_call
    )

    logger.debug(
        'Completed FastQC for %s reads for sample: %s',
        level, sample.name
    )
    return out_str, err_str


def _run_command(
    *,  # Enforce keyword arguments
    command: str,
    logger: logging.Logger
) -> Tuple[str, str]:
    """
    Run a subprocess command and return the output and error.

    Args:
        command (str): The command to run.
        logger (logging.Logger): Logger for recording errors.

    Returns:
        tuple: Standard output and error strings.
    """
    logger.debug('Running command: %s', command)
    out, err = run_subprocess(command=command)
    return out, err


def _log_fastqc_output(
    *,  # Enforce keyword arguments
    err_str: str,
    fastqc_reads: str,
    log_file: str,
    logger: logging.Logger,
    out_str: str,
    sample: CustomBox,
    system_call: str,
) -> None:
    """
    Log the output and errors from the FastQC commands.

    Args:
        system_call (str): FastQC system call.
        fastqc_reads (str): FastQC reads call.
        out_str (str): Standard output string.
        err_str (str): Standard error string.
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
    """
    logger.debug('Logging FastQC output for sample: %s', sample.name)
    write_to_log_file(
        out=system_call,
        err=system_call,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )
    write_to_log_file(
        out=fastqc_reads,
        err=fastqc_reads,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )
    write_to_log_file(
        out=out_str,
        err=err_str,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )


def rename_fastqc_outputs(
    *,  # Enforce keyword arguments
    level: str,
    logger: logging.Logger,
    output_dir: str,
    sample: CustomBox
) -> None:
    """
    Rename the FastQC output files.

    Args:
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        logger (logging.Logger): Logger for recording errors.
        output_dir (str): Output directory.
        sample (CustomBox): Metadata sample object.
    """
    try:
        # Rename the FastQC HTML output file
        shutil.move(
            src=os.path.join(output_dir, 'stdin_fastqc.html'),
            dst=os.path.join(output_dir, f'{sample.name}_fastqc.html')
        )
        # Rename the FastQC ZIP output file
        shutil.move(
            src=os.path.join(output_dir, 'stdin_fastqc.zip'),
            dst=os.path.join(output_dir, f'{sample.name}_fastqc.zip')
        )
        logger.debug(
            'Renamed FastQC output files for %s reads for sample: %s',
            level, sample.name
        )
    except IOError:
        logger.error(
            'Failed to rename FastQC output files for %s reads for sample: %s',
            level, sample.name
        )
