#!/usr/bin/env python3

"""
Script to move FASTQ files for each sample to an appropriately named folder.
"""

# Standard imports
from glob import glob
import logging
import os
import re
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import make_path
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


class FastqMover:
    """
    Class to handle the movement of FASTQ files for each sample to an
    appropriately named folder.
    """

    def __init__(
        self,
        *,  # Enforce keyword arguments
        logger: logging.Logger,
        metadata: List[CustomBox],
        path: str
    ) -> None:
        """
        Initialize the FastqMover class.

        :param logger: Logger for recording information.
        :param metadata: List of metadata sample objects.
        :param path: Path to the directory containing the FASTQ files.
        """
        self.metadata = metadata
        self.path = path
        self.logger = logger
        self.move_fastq()

    def move_fastq(self) -> None:
        """
        Find .fastq files for each sample and move them to an appropriately
        named folder.
        """
        self.logger.info('Moving FASTQ files')
        # Iterate through each sample
        for sample in self.metadata:
            # Retrieve the output directory
            output_dir = os.path.join(self.path, sample.name)

            # Find any fastq files with the sample name
            fastq_files = self.find_fastq_files(
                sample_name=sample.name
            )

            # Only try and move the files if the files exist
            if fastq_files:
                make_path(output_dir)

                # Symlink the fastq files to the directory
                self.symlink_fastq_files(
                    fastq_files=fastq_files,
                    logger=self.logger,
                    output_dir=output_dir
                )

                # Filter out unwanted fastq files
                fastq_files = self.filter_fastq_files(
                    output_dir=output_dir,
                    sample_name=sample.name
                )
            else:
                if output_dir:
                    # Filter out unwanted fastq files
                    fastq_files = self.filter_fastq_files(
                        output_dir=output_dir,
                        sample_name=sample.name
                    )
            sample.general.fastq_files = fastq_files

    def find_fastq_files(self, *, sample_name: str) -> List[str]:
        """
        Find .fastq files for a given sample name.

        This method searches for .fastq files that match the sample name
        using a regular expression. It returns the list of fastq file paths.

        :param sample_name: Name of the sample.
        :return: List of fastq file paths.
        """
        # Compile a regular expression to match all desired patterns
        pattern = re.compile(
            rf'{re.escape(sample_name)}(_.*|\.|.*)\.(fastq|fq)(\.gz)?$'
        )
        # Find all files in the directory
        all_files = glob(os.path.join(self.path, '*'))
        # Filter files using the compiled regular expression
        fastq_files = [
            f for f in all_files if pattern.search(
                os.path.basename(f))]
        return sorted(fastq_files)

    @staticmethod
    def symlink_fastq_files(
        *,
        fastq_files: List[str],
        logger: logging.Logger,
        output_dir: str
    ) -> None:
        """
        Create symlinks for fastq files in the output directory.

        This method creates symbolic links for each fastq file in the specified
        output directory.

        :param fastq_files: List of fastq file paths.
        :param logger: Logger for recording information.
        :param output_dir: Path to the output directory.
        """
        for fastq in fastq_files:
            src = os.path.join('..', os.path.basename(fastq))
            dst = os.path.join(output_dir, os.path.basename(fastq))
            try:
                os.symlink(src, dst)
                logger.info("Created symlink: %s -> %s", src, dst)
            except FileExistsError:
                logger.warning("Symlink already exists: %s", dst)
            except OSError:
                logger.error(
                    "Failed to create symlink: %s -> %s",
                    src,
                    dst,
                    exc_info=True
                )

    @staticmethod
    def filter_fastq_files(*, output_dir: str, sample_name: str) -> List[str]:
        """
        Filter out unwanted fastq files from the output directory.

        This method filters out fastq files that contain unwanted keywords
        such as 'trimmed', 'normalised', 'corrected', 'paired', and 'unpaired'.

        :param output_dir: Path to the output directory.
        :param sample_name: Name of the sample.
        :return: List of filtered fastq file paths.
        """
        # Define unwanted keywords to filter out
        unwanted_keywords = [
            'trimmed', 'normalised', 'corrected', 'paired', 'unpaired'
        ]
        # Find all fastq files in the output directory
        fastq_files = sorted(
            glob(os.path.join(output_dir, f'{sample_name}*.fastq*'))
        )
        # Filter out fastq files containing any of the unwanted keywords
        return [
            fastq for fastq in fastq_files
            if not any(keyword in fastq for keyword in unwanted_keywords)
        ]
