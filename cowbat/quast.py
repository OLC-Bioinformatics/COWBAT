#!/usr/bin/env python3

"""
Assembly evaluation using various tools including Bowtie2, Samtools, Quast,
and Qualimap.
"""

# Standard imports
import logging
import os
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
from typing import List
from io import StringIO

# Third-party imports
from Bio.Sequencing.Applications import SamtoolsIndexCommandline
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file,
    log_str
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def indexing(
    index_queue: Queue,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Use samtools index to index the sorted BAM files.

    Args:
        index_queue (Queue): Queue instance for managing indexing tasks.
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
    logger.info('Indexing sorted BAM files')

    # Guard statements to check preconditions
    if not metadata:
        logger.error("Metadata is empty or not provided.")
        return metadata
    if not isinstance(threads, int) or threads <= 0:
        logger.error("Invalid number of threads specified.")
        return metadata
    if not isinstance(index_queue, Queue):
        logger.error("Index queue is not a Queue instance.")
        return metadata

    # Use ThreadPoolExecutor to manage threads
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit the index method to the executor
        for _ in range(threads):
            executor.submit(index, index_queue)

        # Enqueue tasks
        for sample in metadata:
            if sample.general.best_assembly_file == 'NA':
                continue

            bam_index = SamtoolsIndexCommandline(
                input=sample.quast.sorted_bam
            )
            sample.quast.sorted_bai = sample.quast.sorted_bam + '.bai'
            sample.quast.bam_index = str(bam_index)
            logger.debug(
                "Enqueuing indexing task for sample: %s", sample.name
            )
            index_queue.put((sample, bam_index))

        # Add sentinel values to the queue to signal the threads to exit
        for _ in range(threads):
            index_queue.put(None)

        # Wait for all tasks to be completed
        index_queue.join()

    return metadata


def index(index_queue: Queue) -> None:
    """
    Worker method to process the indexing tasks.

    Args:
        index_queue (Queue): Queue instance for managing indexing tasks.
    """
    while True:
        task = index_queue.get()
        if task is None:
            # Sentinel value encountered, exit the loop
            index_queue.task_done()
            break

        sample, bam_index = task

        # Only make the call if the .bai file doesn't already exist
        if os.path.isfile(sample.quast.sorted_bai):
            index_queue.task_done()
            continue

        # Use StringIO streams to handle output
        stdout, stderr = map(
            StringIO, bam_index(cwd=sample.quast.output_dir)
        )
        if stdout.getvalue() or stderr.getvalue():
            # Set the name of the samtools log file
            samtools_log = os.path.join(
                sample.quast.output_dir,
                'indexing_samtools_bam_index.log'
            )

            # Write the standard error to log
            with open(samtools_log, 'a+', encoding='utf-8') as log:
                log.writelines(
                    log_str(
                        bam_index,
                        stderr.getvalue(),
                        stdout.getvalue()
                    )
                )
        stderr.close()
        index_queue.task_done()


def run_quast(
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    threads: int
) -> List[CustomBox]:
    """
    Run quast on the samples.

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
    logger.info('Running Quast on assemblies')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        sample = prepare_quast_command(
            logger=logger,
            sample=sample,
            threads=threads
        )

        if not quast_report_exists(sample):
            logger.info(
                "Quast report not found for sample %s, running command",
                sample.name
            )
            run_quast_command(
                log_file=log_file,
                logger=logger,
                sample=sample,
            )
        else:
            logger.info(
                "Quast report already exists for sample %s, skipping",
                sample.name
            )

    return metadata


def prepare_quast_command(
    logger: logging.Logger,
    sample: CustomBox,
    threads: int
) -> CustomBox:
    """
    Prepare the Quast command for the sample.

    Args:
        logger (logging.Logger): Logger for recording information.
        sample (CustomBox): Metadata sample object.
        threads (int): Number of threads to use.

    Returns:
        sample (CustomBox): Metadata sample object with quast command
    """
    sample.quast.report = os.path.join(sample.quast.output_dir, 'report.tsv')
    if sample.general.best_assembly_file == "NA":
        logger.warning(
            "Skipping sample %s as best assembly file is 'NA'", sample.name
        )
        return sample

    if len(sample.general.trimmed_corrected_fastq_files) == 2:
        sample.quast.cmd = (
            f'quast.py --pe1 '
            f'{sample.general.trimmed_corrected_fastq_files[0]} '
            f'--pe2 '
            f'{sample.general.trimmed_corrected_fastq_files[1]}'
        )
    else:
        sample.quast.cmd = (
            f'quast.py --single '
            f'{sample.general.trimmed_corrected_fastq_files[0]}'
        )

    sample.quast.cmd += (
        f' --ref-bam {sample.quast.sorted_bam} -t {threads} '
        f'--k-mer-stats --circos --rna-finding '
        f'--conserved-genes-finding -o {sample.quast.output_dir} '
        f'--debug {sample.general.assembly_file} --threads {threads}'
    )
    logger.debug("Constructed quast command: %s", sample.quast.cmd)

    return sample


def quast_report_exists(sample: CustomBox) -> bool:
    """
    Check if the Quast report already exists.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        bool: True if the Quast report exists, False otherwise.
    """
    return os.path.isfile(sample.quast.report)


def run_quast_command(
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Run the Quast command for the sample.

    Args:
        sample (CustomBox): Metadata sample object.
        log_file (str): Path to the log file.

    Raises:
        RuntimeError: If the subprocess command fails.
    """
    try:
        out, err = run_subprocess(sample.quast.cmd)
        logger.debug("Command output: %s", out)
        logger.debug("Command error: %s", err)

        # Write the appropriate information to the log_file
        write_to_log_file(
            out=f'{sample.quast.cmd}\n{out}',
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except RuntimeError as exc:
        logger.error(
            "Failed to run quast command for sample %s: %s", sample.name, exc
        )


def parse_quast_report(
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Parse the quast report, and populate the metadata object with the
    extracted key: value pairs.

    Args:
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata sample objects.

    Raises:
        FileNotFoundError: If the quast report file is not found.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Parsing Quast reports')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        if not os.path.isfile(sample.quast.report):
            logger.warning(
                "Quast report not found for sample %s, skipping",
                sample.name
            )
            continue

        logger.info(
            "Quast report found for sample %s, parsing report",
            sample.name
        )
        sample = _parse_report(
            logger=logger,
            sample=sample
        )

    return metadata


def _parse_report(
    logger: logging.Logger,
    sample: CustomBox
) -> CustomBox:
    """
    Parse the quast report file and extract key-value pairs.

    Args:
        logger (logging.Logger): Logger for recording information
        sample (CustomBox): The sample object to parse the report for.

    Returns:
        CustomBox: Updated metadata sample object.
    """
    try:
        with open(sample.quast.report, 'r', encoding='utf-8') as report:
            for line in report:
                key, value = analyze(line)
                logger.debug(
                    "Extracted key-value pair: %s: %s", key, value
                )
                setattr(sample.quast, key, value)
    except FileNotFoundError as exc:
        logger.error(
            "Quast report file not found for sample %s: %s",
            sample.name, exc
        )

    return sample


def clean_quast(
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> None:
    """
    Remove all the unnecessary temporary files created by quast.

    Args:
        logger (logging.Logger): Logger for recording information
        metadata (List[CustomBox]): List of metadata sample objects.

    Raises:
        FileNotFoundError: If any required file paths are not found.
    """
    logger.info('Cleaning Quast temporary files')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        if not os.path.isdir(sample.quast.output_dir):
            logger.warning(
                "Quast output directory not found for sample %s, skipping",
                sample.name
            )
            continue

        _clean_sample_quast_files(
            logger=logger,
            sample=sample
        )


def _clean_sample_quast_files(
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Clean unnecessary temporary files for a sample.

    Args:
        logger (logging.Logger): Logger for recording information
        sample (CustomBox): The sample object to clean files for.
    """
    for path, _, files in os.walk(sample.quast.output_dir):
        for quast_file in files:
            file_path = os.path.join(path, quast_file)
            _remove_large_files(
                file_path=file_path,
                logger=logger,
                quast_file=quast_file,
                sample=sample
            )


def _remove_large_files(
    file_path: str,
    logger: logging.Logger,
    quast_file: str,
    sample: CustomBox
) -> None:
    """
    Remove large files that do not have .err or .html extensions.

    Args:
        file_path (str): The absolute path of the file.
        logger (logging.Logger): Logger for recording information.
        quast_file (str): The name of the file.
        sample (CustomBox): The sample object to log information for.
    """
    try:
        if (
            os.path.getsize(file_path) > 100000
            and '_sorted.bam' not in quast_file
            and '.err' not in quast_file
            and '.html' not in quast_file
        ):
            logger.info(
                "Removing file %s for sample %s", file_path, sample.name
            )
            os.remove(file_path)
    except FileNotFoundError as exc:
        logger.error(
            "File not found %s for sample %s: %s", file_path, sample.name, exc
        )
    except OSError as exc:
        logger.error(
            "Error removing file %s for sample %s: %s", file_path, sample.name,
            exc
        )


def analyze(line):
    """
    Analyze a line from the report.

    Args:
        line: A line from the report.

    Returns:
        A tuple containing the sanitized key and value.
    """
    # Split the line on tab character
    key, value = line.rstrip().split('\t', 1)

    # Sanitize the key
    key = (
        key.replace(' (%)', '')
        .replace(' ', '_')
        .replace('#', 'num')
        .replace('(>=', 'greater_than')
        .replace(')', '')
        .replace('.', '')
        .replace('\'', '')
    )

    return key, value
