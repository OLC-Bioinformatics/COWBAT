#!/usr/env/python3

"""
Filter contigs from assemblies based on sequence depth and contig length
"""

# Standard imports
from concurrent.futures import as_completed, ThreadPoolExecutor
import logging
import os
import re
import shutil
from typing import List

# Third-party imports
from Bio import SeqIO
from olctools.accessoryFunctions.metadata import CustomBox


def filter_contigs(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    metadata: List[CustomBox],
    sequence_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Filter contigs based on depth and length.

    This method initializes threads and enqueues tasks for filtering
    contigs. It waits for all tasks to be completed.

    Arguments:
    - `logger` should be the logger object.
    - `metadata` should contain the sample data with necessary
        attributes.
    - 'sequence_path' the name and path of the folder with the sequence files
    - `threads` should be set to the number of threads to use.
    """
    logger.info('Filtering contigs')

    futures = []

    # Use ThreadPoolExecutor to manage threads
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Enqueue tasks for each sample
        for sample in metadata:
            if _is_valid_filter_sample(logger=logger, sample=sample):
                logger.debug(
                    "Initializing filter for sample: %s", sample.name
                )
                sample = _initialize_sample_filter(
                    logger=logger,
                    sample=sample
                )
                logger.debug(
                    "Submitting filter task for sample: %s", sample.name
                )
                future = executor.submit(
                    _run_filter,
                    logger=logger,
                    sample=sample,
                    sequence_path=sequence_path
                )
                futures.append(future)

    # Collect results from futures
    updated_metadata = []
    for future in as_completed(futures):
        try:
            sample = future.result()
            logger.debug("Filter task completed for sample: %s", sample.name)
            updated_metadata.append(sample)
        except FileNotFoundError as exc:
            logger.error("File not found during filtering: %s", exc)
        except AttributeError as exc:
            logger.error("Attribute error during filtering: %s", exc)

    return updated_metadata


def _is_valid_filter_sample(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> bool:
    """
    Check if the sample is valid for filtering.

    Args:
        logger: The logger object.
        sample: The sample object to check.

    Returns:
        bool: True if the sample is valid, False otherwise.
    """
    logger.debug(
        "Checking if sample is valid for filtering: %s",
        sample.name)

    # Check if the best assembly file is not 'NA'
    if sample.general.best_assembly_file == 'NA':
        logger.debug(
            "Sample %s is not valid for filtering: best_assembly_file is "
            "'NA'", sample.name
        )
        return False

    # Check if the necessary attributes are present
    if not hasattr(
            sample,
            'general') or not hasattr(
            sample.general,
            'assembly_file'):
        logger.debug(
            "Sample %s is not valid for filtering: missing assembly_file",
            sample.name
        )
        return False
    return True


def _initialize_sample_filter(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Initialize the filter attributes for a sample.

    Args:
        logger: The logger object.
        sample: The sample object to initialize attributes for.
    """
    logger.debug(
        "Initializing filter attributes for sample: %s",
        sample.name
    )
    # Set the name of the unfiltered assembly output file
    sample.general.contigs_file = sample.general.assembly_file

    return sample


def _run_filter(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    sample: CustomBox,
    sequence_path: str
) -> None:
    """
    Run the filter process for a sample.

    Args:
        logger: The logger object.
        sample: The sample object to run the filter on.
        sequence_path: The name and path of the folder containing the
        sequence files
    """
    # Check if the contigs file exists
    if os.path.isfile(sample.general.contigs_file):
        logger.debug("Contigs file exists for sample: %s", sample.name)

        # Filter the contigs based on depth and length
        _filter_contigs(
            logger=logger,
            sample=sample
        )

        # Copy the filtered file to the BestAssemblies folder
        sample = _copy_filtered_file(
            logger=logger,
            sample=sample,
            sequence_path=sequence_path
        )
    else:
        logger.debug(
            "Contigs file does not exist for sample: %s",
            sample.name)

        # Set the best assembly file to 'NA' if contigs file doesn't exist
        sample.general.best_assembly_file = 'NA'

    return sample


def _filter_contigs(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    sample: CustomBox
) -> None:
    """
    Filter contigs based on depth and length for a sample.

    Args:
        logger: The logger object.
        sample: The sample object to filter contigs for.
    """
    logger.debug("Filtering contigs for sample: %s", sample.name)
    pass_depth = []

    # Open the contigs file for reading
    with open(sample.general.contigs_file, 'r', encoding='utf-8') as file:
        # Parse the contigs file in FASTA format
        for record in SeqIO.parse(file, "fasta"):
            # Extract the contig name without '_pilon'
            contig = record.id.split('_pilon')[0]
            # Get the mean coverage and standard deviation
            coverage_mean = float(sample.qualimap.coverage[contig])
            coverage_std = float(sample.qualimap.std_dev[contig])

            # Check if the contig passes the depth and length filters
            if (float(sample.qualimap.coverage[contig]) >
                (coverage_mean - coverage_std * 1.5) and
                    len(record.seq) > 500):

                # Replace 'Contig' in the record ID with the sample name
                new_id = re.sub("Contig", sample.name, record.id)
                record.id = new_id
                # Clear the name and description attributes
                record.name = ''
                record.description = ''
                # Add the record to the list of passing contigs
                pass_depth.append(record)

    # Check if there are any contigs that passed the filters
    if pass_depth:
        logger.debug(
            "Contigs passed the filters for sample: %s",
            sample.name)
        # Set the name of the filtered file
        filtered = sample.general.filtered_file

        # Open the filtered file for writing
        with open(filtered, 'w', encoding='utf-8') as formatted:
            # Write the passing contigs to the filtered file
            SeqIO.write(pass_depth, formatted, 'fasta')


def _copy_filtered_file(
    *,  # Enforce the use of keyword arguments
    logger: logging.Logger,
    sample: CustomBox,
    sequence_path: str
) -> None:
    """
    Copy the filtered file to the BestAssemblies folder.

    Args:
        sample: The sample object to copy the filtered file for.
        sequence_path: The name and path of the folder containing the sample-
            specific sequence files
    """
    logger.debug("Copying filtered file for sample: %s", sample.name)
    if os.path.isfile(sample.general.filtered_file):
        sample.general.best_assemblies_path = os.path.join(
            sequence_path, 'BestAssemblies'
        )
        best_assembly_file = os.path.join(
            sample.general.best_assemblies_path, f'{sample.name}.fasta'
        )
        sample.general.best_assembly_file = best_assembly_file

        if not os.path.isfile(best_assembly_file):
            logger.debug(
                "Copying file to BestAssemblies for sample: %s",
                sample.name
            )
            shutil.copyfile(
                sample.general.filtered_file,
                best_assembly_file)
    else:
        logger.debug(
            "Filtered file does not exist for sample: %s",
            sample.name
        )
        sample.general.best_assembly_file = 'NA'

    return sample
