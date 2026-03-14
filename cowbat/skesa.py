#!/usr/bin/env python3

"""
Run SKESA on FASTQ files.
"""

# Standard imports
import logging
import os
import shutil
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


class Skesa:
    """
    Assemble genomes with SKESA, and ensure that the assembly was successful.
    """

    def __init__(
        self,
        *,  # Enforce keyword arguments
        log_file: str,
        logger: logging.Logger,
        metadata: List[CustomBox],
        report_path: str,
        sequence_path: str,
        threads: int
    ) -> None:
        """
        Initialize the Skesa class with the provided arguments.

        Args:
            log_file (str): Path to the log file.
            logger (logging.Logger): Logger for recording information.
            metadata (List[CustomBox]): List of metadata sample objects.
            report_path (str): Path to the report directory.
            sequence_path (str): Path to the sequence directory.
            threads (int): Number of threads to use.
        """
        self.metadata = metadata
        self.threads = threads
        self.path = sequence_path
        self.logfile = log_file
        self.logger = logger
        self.reportpath = report_path

        # Create the directories as required
        os.makedirs(os.path.join(self.path, 'BestAssemblies'), exist_ok=True)
        os.makedirs(os.path.join(self.path, 'raw_assemblies'), exist_ok=True)
        os.makedirs(self.reportpath, exist_ok=True)
        logger.info('Assembling sequences')

    def main(self) -> List[CustomBox]:
        """
        Assemble genomes and determine whether the assembly was successful.
        """
        self.skesa_assemble()
        self.best_assembly_file()

        return self.metadata

    def skesa_assemble(self) -> None:
        """
        Run SKESA to assemble genomes.
        """
        for sample in self.metadata:
            self.logger.debug("Processing sample: %s", sample.name)
            if not self._initialize_sample_assembly(sample=sample):
                continue

            if sample.commands.assemble and not os.path.isfile(
                sample.general.assembly_file
            ):
                self._run_skesa(sample=sample)

    def _initialize_sample_assembly(self, *, sample: CustomBox) -> bool:
        """
        Initialize the assembly attributes for a sample.

        Args:
            sample: The sample object to initialize attributes for.

        Returns:
            bool: True if initialization is successful, False otherwise.
        """
        try:
            if hasattr(sample.general, 'trimmed_corrected_fastq_files'):
                fastq_files = sample.general.trimmed_corrected_fastq_files
                if not isinstance(fastq_files, list):
                    self.logger.warning(
                        'No FASTQ files found for sample: %s',
                        sample.name
                    )
                    return False

                # Set the output directory
                sample.general.assembly_output = os.path.join(
                    sample.general.output_directory, 'assembly_output'
                )
                os.makedirs(sample.general.assembly_output, exist_ok=True)
                sample.general.assembly_file = os.path.join(
                    sample.general.assembly_output,
                    f'{sample.name}_unfiltered.fasta'
                )
                sample.general.best_assembly_file = os.path.join(
                    sample.general.assembly_output,
                    f'{sample.name}.fasta'
                )

                # Set the forward fastq files
                sample.general.assembly_fastq = fastq_files
                forward = fastq_files[0]
                gz = '.gz' in forward

                # Construct the SKESA command
                if len(fastq_files) == 2:
                    sample.commands.assemble = (
                        f'skesa --fastq {",".join(fastq_files)} '
                        f'--cores {self.threads} --use_paired_ends '
                        f'--vector_percent 1 --contigs_out '
                        f'{sample.general.assembly_file}'
                    )
                else:
                    sample.commands.assemble = (
                        f'skesa --fastq {",".join(fastq_files)} '
                        f'--cores {self.threads} --vector_percent 1 '
                        f'--contigs_out {sample.general.assembly_file}'
                    )

                # Specify that the files are gzipped
                if gz:
                    sample.commands.assemble += ' --gz'
                self.logger.debug(
                    "SKESA command for sample %s: %s",
                    sample.name,
                    sample.commands.assemble)
                return True
            else:
                sample.general.assembly_output = 'NA'
                sample.general.assembly_fastq = 'NA'
                sample.general.best_assembly_file = 'NA'
                return False
        except AttributeError as exc:
            self.logger.error(
                'AttributeError for sample %s: %s', sample.name, str(exc)
            )
            sample.general.assembly_output = 'NA'
            sample.general.assembly_fastq = 'NA'
            sample.general.trimmed_corrected_fastq_files = 'NA'
            sample.general.best_assembly_file = 'NA'
            return False

    def _run_skesa(self, *, sample: CustomBox) -> None:
        """
        Run the SKESA command for a sample.

        Args:
            sample: The sample object to run SKESA on.
        """
        self.logger.debug("Running SKESA for sample: %s", sample.name)
        out, err = run_subprocess(sample.commands.assemble)
        self._log_skesa_output(sample, out, err)

    def _log_skesa_output(self, sample, out, err) -> None:
        """
        Log the output and errors from the SKESA command.

        Args:
            sample: The sample object to log information for.
            out: The standard output from the SKESA command.
            err: The standard error from the SKESA command.
        """
        self.logger.debug("Logging SKESA output for sample: %s", sample.name)
        write_to_log_file(
            out=sample.commands.assemble,
            err=sample.commands.assemble,
            log_file=self.logfile,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=self.logfile,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err
        )

    def best_assembly_file(self) -> None:
        """
        Determine whether the contigs.fasta output file from the assembler is
        present. If not, set the .best_assembly_file attribute to 'NA'.
        """
        for sample in self.metadata:
            self.logger.debug(
                "Checking best assembly file for sample: %s",
                sample.name)
            try:
                self._check_and_copy_assembly_file(sample)
            except AttributeError as e:
                self.logger.error(
                    'AttributeError for sample %s: %s', sample.name, str(e)
                )
                sample.general.assembly_file = 'NA'
                sample.general.best_assembly_file = 'NA'

    def _check_and_copy_assembly_file(self, sample) -> None:
        """
        Check and copy the assembly file for a sample.

        Args:
            sample: The sample object to check and copy the assembly file for.
        """
        # Set the name of the filtered assembly file
        filtered_output_file = os.path.join(
            self.path, 'raw_assemblies', f'{sample.name}.fasta'
        )
        # Set the name of the unfiltered SKESA assembly output file
        if os.path.isfile(sample.general.assembly_file):
            size = os.path.getsize(sample.general.assembly_file)
            # Ensure that the assembly isn't just an empty file
            if size == 0:
                sample.general.best_assembly_file = 'NA'
            else:
                sample.general.best_assembly_file = (
                    sample.general.assembly_file
                )
                shutil.copyfile(
                    sample.general.best_assembly_file,
                    filtered_output_file
                )
        else:
            sample.general.best_assembly_file = 'NA'
        # Add the name and path of the filtered file to the metadata
        sample.general.filtered_file = filtered_output_file
