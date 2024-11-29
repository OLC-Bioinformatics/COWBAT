#!/usr/bin/env python3

"""
A reduced version of the COWBAT pipeline for diagnostic laboratories
"""

# Standard imports
from argparse import ArgumentParser
import inspect
import multiprocessing
import os
import traceback
from typing import (
    Callable,
    List,
    Tuple
)

# Local imports
from cowbat.assemble import assemble
from cowbat.methods import (
    check_programs,
    initialize_logging,
    read_checkpoint,
    sample_metadata,
    tilde_expand,
    write_checkpoint,
)
from cowbat.quality import quality
from cowbat.quality_report import write_quality_report
from cowbat.teacup_version import __version__
from cowbat.taxonomy import taxonomy
from cowbat.typing import typing

__author__ = 'adamkoziol'


class TeacupCOWBAT:
    """
    Runs the minimal COWBAT pipeline on sequence data. Determines checkpoints
    if necessary, stores log files, metadata, and reports.
    """

    def __init__(
        self,
        sequence_path: str,
        database_path: str,
        logging_level: str,
        threads: int = None
    ):
        # Initialise the variables
        self.sequence_path = tilde_expand(path=sequence_path)
        self.database_path = tilde_expand(path=database_path)
        self.report_path = os.path.join(self.sequence_path, 'reports')

        # Set the name of the checkpoint file
        self.checkpoint_file = os.path.join(
            self.sequence_path, '.checkpoint'
        )

        # Set the name of the log files
        self.log_file = os.path.join(self.sequence_path, 'log_file')
        self.error_log_file = os.path.join(self.sequence_path, 'error.log')

        # Initialize logging
        self.logger, self.error_logger = initialize_logging(
            error_log_file=self.error_log_file,
            log_file=self.log_file,
            logging_level=logging_level
        )

        # Print a welcome message
        self.logger.info(
            "Welcome to Teacup COWBAT version %s", __version__
        )

        # Set the list of required command line programs
        programs = [
            'bbduk.sh',
            'bowtie2',
            'centrifuge',
            'confindr',
            'fastqc',
            'mash',
            'metaphlan',
            'multiqc',
            'pilon',
            'prodigal',
            'qualimap',
            'quast.py',
            'samtools',
            'skesa'
        ]

        # Ensure that all the required programs are present in the environment
        missing = check_programs(
            logger=self.logger,
            programs=programs
        )

        # If there are missing programs, exit
        if missing:
            raise SystemExit

        # Determine the last successful step
        self.last_step = read_checkpoint(checkpoint_file=self.checkpoint_file)
        self.logger.info("Starting from step: %s", self.last_step)

        # Define the pipeline steps and their corresponding methods
        self.pipeline_steps: List[Tuple[str, Callable[[], None]]] = [
            ('quality', quality),
            ('assemble', assemble),
            ('taxonomy', taxonomy),
            ('quality_report', write_quality_report),
            ('typing', typing)
        ]

        # Initialise the list of metadata
        self.metadata = []

        # Use the argument for the number of threads to use, or default to the
        # number of cpus in the system
        self.threads = threads if threads else multiprocessing.cpu_count() - 1

    def main(self):
        """
        Run the Teacup COWBAT methods
        """
        self.metadata = sample_metadata(
            error_logger=self.error_logger,
            logger=self.logger,
            metadata=self.metadata,
            sequence_path=self.sequence_path
        )
        self.run_next_step()

    def run_next_step(self) -> None:
        """
        Run the next step in the pipeline based on the last successful step.
        """
        start_index = 0

        # Determine the starting index based on the last successful step
        if self.last_step:
            for i, (step, _) in enumerate(self.pipeline_steps):
                if step == self.last_step:
                    start_index = i + 1
                    break

        # Execute the steps from the determined starting index
        for step, method in self.pipeline_steps[start_index:]:

            self.logger.info('Running %s section of COWBAT', step)

            try:
                # Get the method's parameters
                params = inspect.signature(method).parameters

                # Prepare the arguments to pass to the method
                args = {
                    'analysis_type': step,
                    'error_logger': self.error_logger,
                    'commit': __version__,
                    'log_file': self.log_file,
                    'logger': self.logger,
                    'metadata': self.metadata,
                    'reference_file_path': self.database_path,
                    'report_path': self.report_path,
                    'sequence_path': self.sequence_path,
                    'threads': self.threads
                }

                # Filter the arguments based on the method's parameters
                filtered_args = {k: v for k, v in args.items() if k in params}

                # Call the method with the filtered arguments
                self.metadata = method(**filtered_args)

                # Write the checkpoint
                write_checkpoint(
                    checkpoint_file=self.checkpoint_file,
                    step=step
                )
                self.logger.info("Saved checkpoint at step: %s", step)
            except Exception as exc:
                self.logger.error(
                    "Pipeline terminated due to an error at step '%s': %s",
                    step, exc
                )
                self.logger.error(traceback.format_exc())
                self.error_logger.error(traceback.format_exc())
                raise


def cli():
    """
    Command line argument parser
    """
    # Parser for arguments
    parser = ArgumentParser(
        description='Assemble genomes from Illumina fastq files')
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'%(prog)s commit {__version__}'
    )
    parser.add_argument(
        '-s', '--sequence_path',
        required=True,
        help='Path of the folder containing sequencing reads'
    )
    parser.add_argument(
        '-d', '--database_path',
        required=True,
        help='Path of the folder containing the database'
    )
    parser.add_argument(
        '-t', '--threads',
        help='Number of threads. Default is the number of cores in the system'
    )
    parser.add_argument(
        '-l', '--log',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Set the logging verbosity level'
    )

    # Parse the arguments
    arguments = parser.parse_args()

    # Create the teacup object
    teacup = TeacupCOWBAT(
        sequence_path=arguments.sequence_path,
        database_path=arguments.database_path,
        logging_level=arguments.log,
        threads=arguments.threads
    )

    # Run the Teacup COWBAT pipeline
    teacup.main()


if __name__ == '__main__':
    cli()
