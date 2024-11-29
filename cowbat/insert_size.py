#!/usr/bin/env python3

"""
Estimate library insert size
"""

# Standard imports
import logging
import os
from typing import List

from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def extract_insert_size(
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Parse the bwa index log information to extract insert size estimations.

    Args:
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.

    Raises:
        FileNotFoundError: If the reads stats file is not found.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Calculating insert size')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        sample.quast.reads_stats_file = os.path.join(
            sample.quast.output_dir, 'reads_stats', 'reads_stats.err'
        )

        if os.path.isfile(sample.quast.reads_stats_file):
            logger.info(
                "Reads stats file found for sample %s, parsing file",
                sample.name
            )

            sample = initialize_sample_attributes(sample)
            sample = parse_reads_stats_file(sample)
        else:
            logger.warning(
                "Reads stats file not found for sample %s, skipping",
                sample.name
            )
            sample.quast.insert_mean = 'ND'
            sample.quast.insert_std = 'ND'

    return metadata


def initialize_sample_attributes(sample: CustomBox) -> CustomBox:
    """
    Initialize attributes for the insert size estimation.

    Args:
        sample (CustomBox): The sample object to initialize attributes for.

    Returns:
        sample (CustomBox): The updated sample object
    """
    sample.quast.total_reads = 0
    sample.quast.insert_mean = []
    sample.quast.insert_std = []
    sample.quast.read_blocks = []

    return sample


def parse_reads_stats_file(sample: CustomBox) -> CustomBox:
    """
    Parse the reads stats file to extract insert size estimations.

    Args:
        sample (CustomBox): The sample object to populate with extracted data.

    Returns:
        sample (CustomBox): The updated sample
    """
    current_reads = 0

    # Create a variable to store the long attribute name
    stats_file = sample.quast.reads_stats_file

    with open(stats_file, 'r', encoding='utf-8') as read_stats:
        for line in read_stats:

            # Create a variable to store the long FR section line for
            # easier parsing
            size = 'analyzing insert size distribution for orientation FR'

            # BWA estimates the insert size distribution per 256*1024
            # read pairs. Extract the number of reads present in the
            # current block being processed e.g. # candidate unique pairs
            # for (FF, FR, RF, RR): (46,226102, 14, 28)
            if '# candidate unique pairs for' in line:
                current_reads = int(
                    line.rstrip().replace(',', '').replace('(', '')
                    .replace(')', '').split()[-3]
                )
                sample.quast.total_reads += current_reads

            # Continue parsing to find the FR section of the current block
            elif size in line:
                sample = extract_insert_size_distribution(
                    read_stats, sample, current_reads
                )

    return sample


def extract_insert_size_distribution(
    read_stats, sample: CustomBox, current_reads: int
) -> CustomBox:
    """
    Extract the mean and standard deviation of the insert size for a block.

    Args:
        read_stats: The file object for the reads stats file.
        sample (CustomBox): The sample object to populate with extracted data.
        current_reads (int): The number of reads in the current block.

    Returns:
        sample (CustomBox): The updated sample
    """
    for sub_line in read_stats:

        # Extract the mean and standard deviation of the insert size for
        # this block [M::mem_pestat] mean and std.dev: (487.88, 246.14)
        if '[M::mem_pestat] mean and std.dev:' in sub_line:
            split_line = (
                sub_line.rstrip().replace(',', '').replace('(', '')
                .replace(')', '').split()
            )
            mean = float(split_line[-2])
            std = float(split_line[-1])
            sample.quast.insert_mean.append(mean)
            sample.quast.insert_std.append(std)
            sample.quast.read_blocks.append(current_reads)
            break

    return sample


def calculate_weighted_insert_size(
    logger: logging.Logger,
    metadata: List[CustomBox]
) -> List[CustomBox]:
    """
    Calculate the weighted mean and standard deviation of the insert size
    from the extracted bwa blocks.

    Args:
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logger.info('Calculating weighted insert size')

    for sample in metadata:
        logger.debug("Processing sample: %s", sample.name)

        # Initialize attributes to store the calculated weighted mean and
        # standard deviation of insert sizes
        sample.quast.mean_insert = float()
        sample.quast.std_insert = float()

        # Guard statement to check if insert_mean is 'ND'
        if sample.quast.insert_mean == 'ND':
            logger.warning(
                "Insert mean is 'ND' for sample %s, skipping",
                sample.name
            )
            continue

        # Guard statement to check if read_blocks is empty
        if not sample.quast.read_blocks:
            logger.warning(
                "No read blocks found for sample %s, skipping",
                sample.name
            )
            continue

        # Iterate through all the read blocks present in the sample
        for i, read_block in enumerate(sample.quast.read_blocks):

            # Calculate the weight of the current block by dividing it
            # (current number of reads) by the total number of reads in
            # the sample
            weight = read_block / sample.quast.total_reads

            # Multiply the mean for this block to obtain the weighted
            # mean, and add it to the total mean
            sample.quast.mean_insert += (
                sample.quast.insert_mean[i] * weight
            )

            # Same calculation, but for standard deviation
            sample.quast.std_insert += (
                sample.quast.insert_std[i] * weight
            )

        # Set the attributes to floats with two decimal places
        sample.quast.mean_insert = float(
            f'{sample.quast.mean_insert:.2f}'
        )
        sample.quast.std_insert = float(
            f'{sample.quast.std_insert:.2f}'
        )

    return metadata
