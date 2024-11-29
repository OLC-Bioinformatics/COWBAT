#!/usr/bin/env python3
"""
Runs Prodigal for gene prediction on sequence data.
"""

# Standard imports
from concurrent.futures import ThreadPoolExecutor
import logging
import os
from queue import Queue
from threading import Lock
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'

# Initialise a thread lock for thread-safe logging
thread_lock = Lock()


class Prodigal:
    """
    Class to perform gene predictions using Prodigal.
    """

    def __init__(
        self,
        log_file: str,
        logger: logging.Logger,
        metadata: List[CustomBox]
    ):
        """
        Initialize the Prodigal class.

        Args:
            log_file (str): Path to the log file.
            logger (logging.Logger): Logger for recording information.
            metadata (List[CustomBox]): List of metadata sample objects.

        Preconditions:

        """
        # Extract necessary attributes from the input object
        self.metadata = metadata
        self.log_file = log_file
        self.logger = logger

        self.logger.debug("Initializing Prodigal class")

        # Initialize queues for threading
        self.predict_queue = Queue()
        self.parse_queue = Queue()

    def main(self) -> List:
        """
        Run the prediction and parsing processes
        """
        # Start the prediction and parsing processes
        self.predict_threads()
        self.prodigal_parse()

        # Return the updated metadata
        return self.metadata

    def predict_threads(self) -> None:
        """
        Create threads for gene prediction.

        Preconditions:
        - self.metadata must be a list of sample objects.
        """
        self.logger.info('Performing gene predictions')
        self.logger.debug("Creating threads for gene prediction")

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            for _ in range(os.cpu_count()):
                executor.submit(self.predict)

            # Enqueue tasks
            for sample in self.metadata:
                # Initialize the Prodigal attribute for each sample
                sample.prodigal = CustomBox()
                if sample.general.best_assembly_file == 'NA':
                    continue

                self.logger.debug(
                    "Adding sample to predict queue: %s", sample.name
                )
                # Add the sample to the prediction queue
                self.predict_queue.put(sample)

            # Add sentinel values to the queue to signal the threads to exit
            for _ in range(os.cpu_count()):
                self.predict_queue.put(None)

            # Wait for all tasks in the queue to be completed
            self.predict_queue.join()

    def predict(self) -> None:
        """
        Perform gene prediction on a sample.

        This method runs in a loop, processing tasks from the predict_queue.
        It performs gene prediction and logs any errors.
        """
        while True:
            sample = self.predict_queue.get()
            if sample is None:
                # Sentinel value encountered, exit the loop
                self.predict_queue.task_done()
                break

            self.logger.debug(
                "Starting gene prediction for sample: %s", sample.name
            )

            # Populate attributes for Prodigal results
            self.populate_prodigal_attributes(sample)

            # Create the folder to store the reports
            self.create_report_directory(sample)

            # Check if the report already exists and is not empty
            if self.is_report_needed(sample):
                # Run the Prodigal command
                self.run_prodigal_command(sample)

            self.logger.debug(
                "Completed gene prediction for sample: %s", sample.name
            )

            # Mark the task as done
            self.predict_queue.task_done()

    def populate_prodigal_attributes(self, sample) -> None:
        """
        Populate attributes for Prodigal results.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        sample.prodigal.report_dir = os.path.join(
            sample.general.output_directory, 'prodigal'
        )
        sample.prodigal.results_file = os.path.join(
            sample.prodigal.report_dir,
            f'{sample.name}_prodigalresults.sco'
        )
        sample.prodigal.results = sample.prodigal.results_file

        # Create a variable to store the genes.fa path
        fasta_path = os.path.join(
            sample.prodigal.report_dir, f"{sample.name}_genes.fa"
        )

        sample.commands.prodigal = (
            f'prodigal -i {sample.general.best_assembly_file} '
            f'-o {sample.prodigal.results_file} -f sco -d '
            f'{fasta_path}'
        )

        self.logger.debug(
            "Prodigal command for sample %s: %s",
            sample.name, sample.commands.prodigal
        )

    def create_report_directory(self, sample) -> None:
        """
        Create the folder to store the reports.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        os.makedirs(sample.prodigal.report_dir, exist_ok=True)

    def is_report_needed(self, sample) -> bool:
        """
        Check if the report already exists and is not empty.

        Args:
            sample: A sample object containing metadata and file paths.

        Returns:
            bool: True if the report is needed, False otherwise.
        """
        return not os.path.isfile(sample.prodigal.results_file) or \
            os.stat(sample.prodigal.results_file).st_size == 0

    def run_prodigal_command(self, sample) -> None:
        """
        Run the Prodigal command and log the output.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        out, err = run_subprocess(sample.commands.prodigal)
        with thread_lock:
            # Log the command and its output
            write_to_log_file(
                out=sample.commands.prodigal,
                err=sample.commands.prodigal,
                log_file=self.log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=self.log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

    def prodigal_parse(self) -> None:
        """
        Parse the results of gene predictions.

        Preconditions:
        - self.metadata must be a list of sample objects.
        """
        self.logger.info('Parsing gene predictions')

        for sample in self.metadata:
            # Initialize Prodigal attributes for the sample
            self.initialize_prodigal_attributes(sample)

            if sample.general.best_assembly_file != 'NA':
                self.logger.debug(
                    "Parsing Prodigal results for sample: %s", sample.name
                )
                # Parse the Prodigal results for the sample
                self.parse_prodigal_results(sample)

    def initialize_prodigal_attributes(self, sample) -> None:
        """
        Initialize Prodigal attributes for a sample.

        Args:
            sample: A sample object containing metadata and file paths.

        Preconditions:
        - sample must have a prodigal attribute.
        """

        # Initialize counters for predicted genes
        sample.prodigal.predicted_genes_total = 0
        sample.prodigal.predicted_genes_over_3000bp = 0
        sample.prodigal.predicted_genes_over_1000bp = 0
        sample.prodigal.predicted_genes_over_500bp = 0
        sample.prodigal.predicted_genes_under_500bp = 0

    def parse_prodigal_results(self, sample) -> None:
        """
        Parse the Prodigal results for a sample.

        Args:
            sample: A sample object containing metadata and file paths.

        Preconditions:
        - sample.prodigal.results must be a valid file path.
        """
        self.logger.debug(
            "Opening Prodigal results file for sample: %s", sample.name
        )

        # Open the Prodigal results file for reading
        with open(sample.prodigal.results, 'r', encoding='utf-8') as results:
            for line in results:
                if line.startswith('>'):
                    # Extract the start and end positions of the gene
                    start, end = map(int, line.split('_')[1:3])
                    length = abs(start - end)
                    sample.prodigal.predicted_genes_total += 1

                    # Categorize the gene based on its length
                    if length > 3000:
                        sample.prodigal.predicted_genes_over_3000bp += 1
                    elif length > 1000:
                        sample.prodigal.predicted_genes_over_1000bp += 1
                    elif length > 500:
                        sample.prodigal.predicted_genes_over_500bp += 1
                    else:
                        sample.prodigal.predicted_genes_under_500bp += 1
