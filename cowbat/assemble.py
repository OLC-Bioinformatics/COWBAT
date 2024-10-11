#!/usr/env/python

"""
Collection of functions to assemble bacterial genomes
"""

# Third-party imports
from genemethods.assemblypipeline.assembly_evaluation import AssemblyEvaluation
from genemethods.assemblypipeline.prodigal import Prodigal
from genemethods.assemblypipeline.skesa import Skesa
from olctools.accessoryFunctions.accessoryFunctions import (
    write_metadata_to_file
)

__author__ = 'adamkoziol'


def assemble(
    error_logger,
    log_file,
    metadata,
    report_path,
    sequence_path,
    threads
):
    """
    Assemble genomes and perform some basic quality analyses
    """
    # Assemble genomes
    assembly = Skesa(
        log_file=log_file,
        metadata=metadata,
        report_path=report_path,
        sequence_path=sequence_path,
        threads=threads
    )
    metadata = assembly.main()

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        metadata=metadata
    )

    # Calculate assembly metrics on raw assemblies
    qual = AssemblyEvaluation(
        log_file=log_file,
        metadata=metadata,
        sequence_path=sequence_path,
        threads=threads
    )
    metadata = qual.main()

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        metadata=metadata
    )
    # ORF detection
    prod = Prodigal(
        log_file=log_file,
        metadata=metadata
    )
    metadata = prod.main()

    # Write the metadata to file
    write_metadata_to_file(
        error_logger=error_logger,
        metadata=metadata
    )
    return metadata
