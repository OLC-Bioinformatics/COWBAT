#!/usr/bin/env python3

"""
Assembly evaluation using various tools including Bowtie2, Samtools, Quast,
and Qualimap.
"""

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import os
import multiprocessing
from queue import Queue
from typing import Any, List

# Third-party imports
from cowbat.bowtie2 import (
    bowtie_build,
    bowtie_run
)
from cowbat.filter import filter_contigs
from cowbat.insert_size import (
    calculate_weighted_insert_size,
    extract_insert_size
)
from cowbat.pilon import pilon
from cowbat.qualimap import parse_qualimap_report, qualimapper
from cowbat.quast import (
    clean_quast,
    indexing,
    parse_quast_report,
    run_quast
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


class AssemblyEvaluation:
    """
    Class to evaluate genome assemblies using various tools.
    """

    def __init__(
        self,
        log_file: str,
        logger: logging.Logger,
        metadata: List[Any],
        sequence_path: str,
        threads: int
    ) -> None:
        """
        Initialize the AssemblyEvaluation class with the provided input object.

        Args:
            log_file (str): Path to the log file.
            logger (logging.Logger): Logger for recording information.
            metadata (List[Any]): List of metadata sample objects.
            sequence_path (str): Path to the sequence directory.
            threads (int): Number of threads to use.
        """
        self.metadata = metadata
        self.threads = threads
        self.log_file = log_file
        self.logger = logger
        self.sequence_path = sequence_path

        # Initialize queues
        self.index_queue = Queue(maxsize=self.threads)
        self.filter_queue = Queue(maxsize=self.threads)

    def main(self) -> None:
        """
        Run the methods in the correct order.
        """

        self.logger.info('Running bowtie2 for assembly evaluation')
        self.metadata = bowtie_build(
            log_file=self.log_file,
            logger=self.logger,
            metadata=self.metadata
        )
        self.metadata = bowtie_run(
            log_file=self.log_file,
            logger=self.logger,
            metadata=self.metadata,
            threads=self.threads
        )
        self.metadata = indexing(
            index_queue=self.index_queue,
            logger=self.logger,
            metadata=self.metadata,
            threads=self.threads
        )

        self.logger.info('Running quast and qualimap for assembly evaluation')
        self.metadata = run_quast(
            log_file=self.log_file,
            logger=self.logger,
            metadata=self.metadata,
            threads=self.threads
        )
        self.metadata = parse_quast_report(
            logger=self.logger,
            metadata=self.metadata,
        )
        self.metadata = qualimapper(
            log_file=self.log_file,
            logger=self.logger,
            metadata=self.metadata,
            threads=self.threads
        )
        self.metadata = parse_qualimap_report(
            logger=self.logger,
            metadata=self.metadata
        )
        clean_quast(
            logger=self.logger,
            metadata=self.metadata
        )
        self.metadata = extract_insert_size(
            logger=self.logger,
            metadata=self.metadata
        )
        self.metadata = calculate_weighted_insert_size(
            logger=self.logger,
            metadata=self.metadata
        )
        self.metadata = pilon(
            log_file=self.log_file,
            logger=self.logger,
            metadata=self.metadata,
            threads=self.threads
        )
        self.metadata = filter_contigs(
            logger=self.logger,
            metadata=self.metadata,
            sequence_path=self.sequence_path,
            threads=self.threads
        )
        self.clear()

        return self.metadata

    def clear(self):
        """
        Clear out large attributes from the metadata objects
        """
        for sample in self.metadata:
            try:
                delattr(sample.qualimap, 'bases')
                delattr(sample.qualimap, 'coverage')
                delattr(sample.qualimap, 'length')
                delattr(sample.qualimap, 'stddev')
            except AttributeError:
                pass


class Parser:

    """
    Parser class to allow running the assembly evaluation from the command line
    """
    def __init__(self):
        """
        Initialize the parser and set up the arguments.
        """
        parser = ArgumentParser(
            description='Calculates coverage depth by mapping FASTQ reads '
                        'against assemblies'
        )
        parser.add_argument(
            '-p', '--path',
            default=os.getcwd(),
            help='Specify the path of the folder that either contains the '
            'files of interest, or will be used to store the outputs'
        )
        parser.add_argument(
            '-a', '--assemblies',
            help='Path to a folder of assemblies. If not provided, the script '
            'will look for .fa or .fasta files in the path'
        )
        parser.add_argument(
            '-f', '--fastq',
            help='Path to a folder of fastq files. If not provided, the '
            'script will look for fastq or .fastq.gz files in the path'
        )
        parser.add_argument(
            '-t', '--threads',
            type=int,
            default=multiprocessing.cpu_count(),
            help='Number of threads. Default is the number of cores in the '
            'system'
        )

        # Get the arguments into an object
        args = parser.parse_args()

        # Define variables from the arguments
        self.sequence_path = os.path.join(args.path, '')
        self.assembly_path = os.path.join(
            args.assemblies, '') if args.assemblies else self.sequence_path
        self.fastq_path = os.path.join(
            args.fastq, '') if args.fastq else self.sequence_path
        self.threads = args.threads

        # Initialize variables
        self.strains = []
        self.metadata = []
        self.log_file = os.path.join(self.sequence_path, 'log_file.txt')

        # Associate the assemblies and fastq files in a metadata object
        self.associate()

        # Evaluate the assemblies
        AssemblyEvaluation(
            log_file=self.log_file,
            logger=logging.getLogger(),
            metadata=self.metadata,
            sequence_path=self.sequence_path,
            threads=self.threads
        )

    def associate(self):
        """
        Associate assemblies and fastq files in a metadata object.
        """
        # Get the sequences in the sequences folder into a list
        self.strains = [
            fasta
            for fasta in sorted(
                glob(os.path.join(self.assembly_path, '*.fa*')))
            if '.fastq' not in fasta]

        for strain in self.strains:
            # Extract the name of the strain from the path and file extension
            strain_name = os.path.split(strain)[1].split('.')[0]

            # Find the corresponding fastq files for each strain
            fastq_files = sorted(
                glob(
                    os.path.join(
                        self.fastq_path,
                        f'{strain_name}*fastq*')))

            # Ensure that fastq files are present for each assembly
            assert fastq_files, f'Cannot find fastq files for strain {
                strain_name} '

            # Create the metadata object
            metadata = CustomBox()
            metadata.name = strain_name
            metadata.general = CustomBox()
            metadata.general.best_assemblies_path = self.assembly_path
            metadata.general.trimmed_fastq_files = fastq_files
            metadata.general.output_directory = os.path.join(
                self.sequence_path, strain_name)
            metadata.mapping = CustomBox()

            # Append the metadata for each sample to the list of samples
            self.metadata.append(metadata)


if __name__ == '__main__':
    Parser()
