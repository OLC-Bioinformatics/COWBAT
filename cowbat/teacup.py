#!/usr/bin/env python3

"""
A reduced version of the COWBAT pipeline for diagnostic laboratories
"""

# Standard imports
from argparse import ArgumentParser
import logging
import os
import traceback
from typing import (
    Callable,
    List,
    Tuple
)

# Third-party imports
from genemethods.assemblypipeline.fastqmover import FastqMover

# Local imports
from cowbat.methods import (
    initialize_logging,
    read_checkpoint,
    tilde_expand,
    write_checkpoint,
    write_metadata_to_file
)
from cowbat.set_run_metadata import (
    basic,
    determine_and_parse_sample_sheet,
    process_sample
)
from cowbat.teacup_version import __version__

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
        logging_level: str
    ):
        # Initialise the variables
        self.sequence_path = tilde_expand(path=sequence_path)
        self.database_path = tilde_expand(path=database_path)

        # Set the name of the checkpoint file
        self.checkpoint_file = os.path.join(
            self.sequence_path, '.checkpoint'
        )

        # Set the name of the log files
        self.log_file = os.path.join(self.sequence_path, 'logfile')
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

        # Determine the last successful step
        self.last_step = read_checkpoint(self.checkpoint_file)
        self.logger.info("Starting from step: %s", self.last_step)

        # Define the pipeline steps and their corresponding methods
        self.pipeline_steps: List[Tuple[str, Callable[[], None]]] = [
            ('create_quality_object', self.create_quality_object),
            # ('quality', self.quality),
            # ('assemble', self.assemble),
            # ('agnostic_typing', self.agnostic_typing),
            # ('typing', self.typing)
        ]

        # Initialise the list of metadata
        self.metadata = []

    def main(self):
        """
        Run the Teacup COWBAT methods
        """
        self.sample_metadata()
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
            try:
                method()
                write_checkpoint(self.checkpoint_file, step)
                logging.info("Saved checkpoint at step: %s", step)
            except Exception as exc:
                logging.error(
                    "Pipeline terminated due to an error at step '%s': %s",
                    step, exc
                )
                logging.error(traceback.format_exc())
                self.error_logger(traceback.format_exc())
                raise

    def sample_metadata(self):
        """
        Helper method to prep metadata from sample sheet (if available)
        """
        # Define the sample sheet
        sample_sheet = os.path.join(
            self.sequence_path,
            'SampleSheet.csv'
        )

        # Process the samples appropriate depending on whether a sample sheet
        # was provided
        if os.path.isfile(sample_sheet):

            # Extract the necessary information from the sample sheet
            data = determine_and_parse_sample_sheet(
                file_path=sample_sheet
            )
            # Create the metadata for each sample
            self.metadata = process_sample(
                commit=__version__,
                data=data,
                path=self.sequence_path
            )
        else:
            # If the sample sheet is missing, perform basic assembly
            self.metadata = basic(
                commit=__version__,
                sequence_path=self.sequence_path
            )

        # Move/link the FASTQ files to strain-specific working directories
        FastqMover(
            metadata=self.metadata,
            path=self.sequence_path
        )

        # Write the metadata to file
        write_metadata_to_file(
            error_logger=self.error_logger,
            metadata=self.metadata
        )

    def create_quality_object(self):
        """
        Run quality checking on the samples
        """
        print(self.metadata)
        quit()


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
        logging_level=arguments.log
    )

    # Run the Teacup COWBAT pipeline
    teacup.main()


if __name__ == '__main__':
    cli()
