#!/usr/env/python3

"""
Run GTDB-Tk on the assemblies
"""

# Standard imports
import csv
import logging
from glob import glob
import os
from typing import List


# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    CustomBox,
    relative_symlink,
    run_subprocess,
    write_to_log_file
)


def gtdbtk(
    *,  # Enforce keyword arguments
    assembly_path: str,
    log_file: str,
    logger: logging.Logger,
    reference_file_path: str,
    report_path: str,
    threads: int
) -> str:
    """
    Process the assemblies with GTDB-Tk

    Args:
        assembly_path (str): Path to the assemblies.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.
        threads (int): Number of threads to use for processing.

    Returns:
        str: Path to the GTDB-Tk report.
    Raises:
        RuntimeError: If there is an issue with the GTDB-Tk analyses.
    """

    logger.info('Running GTDB-Tk analyses')

    # Symlink the assemblies to the genomes folder
    _link_assemblies(assembly_path=assembly_path)

    # Create the GTDB-Tk command
    command = _create_gtdbtk_command(
        assembly_path=assembly_path,
        reference_file_path=reference_file_path,
        report_path=report_path,
        threads=threads
    )

    logger.debug('Running GTDB-Tk command: %s', command)

    # See if the GTDB-Tk report already exists
    report = _locate_gtdbtk_report(report_path=report_path)
    if report:
        logger.info('GTDB-Tk report already exists. Skipping GTDB-Tk')
        return report

    # Run GTDB-Tk
    out, err = run_subprocess(
        command=command
    )

    # Write the command to file
    write_to_log_file(
        out=command,
        err=command,
        log_file=log_file
    )

    # Write stdout and stderr to the log file
    write_to_log_file(
        out=out,
        err=err,
        log_file=log_file
    )

    # Locate the GTDB-Tk report
    report = _locate_gtdbtk_report(report_path=report_path)

    return report


def _link_assemblies(
    *,  # Enforce keyword arguments
    assembly_path: str
) -> None:
    """
    Create the genomes subfolder in the BestAssemblies folder, and link the
    assemblies to this folder

    Args:
        assembly_path (str): Path to the assemblies.
    """

    # Create the genomes subfolder
    genomes_path = os.path.join(assembly_path, 'genomes')

    # Create the folder
    os.makedirs(genomes_path, exist_ok=True)

    # Use glob to find all the assemblies
    assemblies = glob(os.path.join(assembly_path, '*.fasta'))

    # Relative symlink the assemblies to the genomes folder.
    for assembly in assemblies:
        relative_symlink(
            src_file=assembly,
            output_dir=genomes_path
        )


def _locate_gtdbtk_report(
    *,  # Enforce keyword arguments
    report_path: str
) -> str:
    """
    Locate the GTDB-Tk report

    Args:
        report_path (str): Path to the directory containing the GTDB-Tk report.

    Returns:
        str: Path to the GTDB-Tk report.
    """

    # Use glob to locate the GTDB-Tk report
    try:
        report = glob(os.path.join(report_path, 'gtdbtk.*summary.tsv'))[0]
    except IndexError:
        report = None

    return report


def _create_gtdbtk_command(
    *,  # Enforce keyword arguments
    assembly_path: str,
    reference_file_path: str,
    report_path: str,
    threads: int
) -> str:
    """
    Create the GTDB-Tk command

    Args:
        assembly_path (str): Path to the assemblies.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.
        threads (int): Number of threads to use for processing.

    Returns:
        str: The GTDB-Tk command.
    """

    # Set the reference database path
    gtdbtk_path = os.path.join(
        reference_file_path, 'gtdbtk', 'release_data'
    )

    # Set the path to the MASH database
    mash_db = os.path.join(reference_file_path, 'mash')

    # Create the GTDB-Tk command
    command = (
        f'docker run --rm -v {assembly_path}:/data '
        f'-v {gtdbtk_path}:/refdata '
        f'-v {report_path}:/reports '
        f'-v {mash_db}:/mash '
        'olcbioinformatics-gtdbtk classify_wf --genome_dir /data/genomes '
        f'--out_dir /reports --mash_db /mash -x fasta --cpus {threads}'
    )

    return command


def parse_gtbdtk_output(
    *,  # Enforce keyword arguments
    metadata: List[CustomBox],
    report: str
) -> List[CustomBox]:
    """
    Parse the GTDB-Tk output

    Args:
        report (str): Path to the GTDB-Tk report.
    """

    # If the report does not exist, raise an error
    if not report:
        raise RuntimeError('GTDB-Tk report not found')

    # Read the report with csv.DictReader
    with open(report, 'r', encoding='utf-8') as report_fh:
        reader = csv.DictReader(report_fh, delimiter='\t')
        for row in reader:
            # Find the sample in the metadata
            for sample in metadata:
                if row['user_genome'] == sample.name:
                    # Create the gtdbtk attribute on the sample (if necessary)
                    if not sample.key_exists('gtbtk'):
                        sample.gtdbtk = CustomBox()

                    # # Add the GTDB-Tk data to the metadata
                    sample.gtdbtk.classification = row['classification']
                    sample.gtdbtk.closest_genome_reference = \
                        row['closest_genome_reference']
                    sample.gtdbtk.closest_genome_reference_radius = \
                        row['closest_genome_reference_radius']
                    sample.gtdbtk.closest_genome_taxonomy = \
                        row['closest_genome_taxonomy']

                    # Extract the taxonomic information
                    sample = extract_taxonomic_info(
                        classification_str=row['classification'],
                        sample=sample
                    )

    # Return the metadata
    return metadata


def extract_taxonomic_info(
    *,  # Enforce keyword arguments
    classification_str: str,
    sample: CustomBox
) -> CustomBox:
    """
    Extract taxonomic information from the classification string

    e.g. d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;.....

    Args:
        classification_str (str): The classification string.
        sample (CustomBox): The sample object.

    Returns:
        CustomBox: The updated sample object.
    """
    # Split the classification string by semicolons
    taxonomic_levels = classification_str.split(';')

    # Initialize a dictionary to store taxonomic information
    taxonomic_info = {
        'd': 'ND',
        'p': 'ND',
        'c': 'ND',
        'o': 'ND',
        'f': 'ND',
        'g': 'ND',
        's': 'ND'
    }

    # Iterate over each taxonomic level
    for level in taxonomic_levels:
        # Split each level by double underscores
        parts = level.split('__')

        # Ensure there are two parts
        if len(parts) != 2:
            continue

        # Extract the prefix and name
        prefix, name = parts

        # Ensure the prefix is in the taxonomic_info dictionary
        if prefix in taxonomic_info:
            # Update the taxonomic_info dictionary
            taxonomic_info[prefix] = name

    # Set the attributes on the sample object
    sample.gtdbtk.domain = taxonomic_info['d']
    sample.gtdbtk.phylum = taxonomic_info['p']
    sample.gtdbtk.class_ = taxonomic_info['c']
    sample.gtdbtk.order = taxonomic_info['o']
    sample.gtdbtk.family = taxonomic_info['f']
    sample.gtdbtk.genus = taxonomic_info['g']
    sample.general.reference_genus = taxonomic_info['g']
    sample.gtdbtk.species = taxonomic_info['s']

    # Return the updated sample
    return sample
