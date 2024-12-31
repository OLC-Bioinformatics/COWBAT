#!/usr/env/python3

"""
Perform rMLST analyses on raw reads using KMA
"""

# Standard imports
import logging
import os
import sqlite3
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import run_subprocess
from olctools.accessoryFunctions.metadata import CustomBox
import pandas as pd

# Local imports
from cowbat.methods import write_to_log_files

__author__ = 'adamkoziol'


def rmlst(
    *,  # Enforce keyword arguments
    error_logger: logging.Logger,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    report_path: str
) -> List[CustomBox]:
    """
    Perform rMLST analyses on raw reads using KMA

    Args:
        error_logger (logging.Logger): Logger for recording errors.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.

    Returns:
        List[CustomBox]: Updated metadata after all processing steps.

    Raises:
        IOError: If there is an issue with file operations.
        RuntimeError: If there is an issue with the assembly or quality
            analyses.
    """
    # Ensure that the SQLite reference database exists
    rmlst_database = _sqlite_database(
        logger=logger,
        reference_file_path=reference_file_path
    )

    # Print the reference genus for each sample
    for sample in metadata:
        # Create the necessary attributes for the metadata objects
        sample = _set_metadata_attributes(
            reference_file_path=reference_file_path,
            sample=sample
        )

        # Run the KMA command if the report does not exist
        if os.path.isfile(sample.rmlst.report):
            continue

        # Log the command
        logger.debug(sample.rmlst.kma_command)

        # Run the KMA command
        out, err = run_subprocess(
            command=sample.commands.rmlst
        )

        # Write the command and the outputs to the log files
        write_to_log_files(
            command=sample.commands.rmlst,
            err=err,
            log_file=log_file,
            logger=logger,
            out=out,
            program='rmlst',
            sample=sample
        )

    # Load the SQLite database
    conn = sqlite3.connect(rmlst_database)

    # Parse the KMA results
    for sample in metadata:
        # Parse the KMA results
        sample = _parse_kma_results(
            logger=logger,
            sample=sample
        )
        # Get matching profiles
        sample = _get_rmlst_profiles(
            conn=conn,
            logger=logger,
            sample=sample
        )

    # Close database connection
    conn.close()

    return metadata


def _sqlite_database(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    reference_file_path: str
) -> str:
    """
    Ensure that the SQLite reference database exists

    Args:
        logger (logging.Logger): Logger object
        reference_file_path (str): Path to the reference database

    Returns:
        str: Path to the SQLite reference database
    """
    # Set the path to the SQLite database
    sqlite_database = os.path.join(
        reference_file_path,
        'rmlst',
        'rmlst.sqlite'
    )

    # Ensure that the SQLite database exists
    if not os.path.isfile(sqlite_database):
        # Create the SQLite database
        profile_path = os.path.join(
            reference_file_path,
            'rmlst',
            'profiles.txt'
        )
        # Create the SQLite database
        _create_rmlst_database(
            database_path=sqlite_database,
            logger=logger,
            profile_path=profile_path
        )
    return sqlite_database


def _create_rmlst_database(
    *,  # Enforce keyword arguments
    database_path: str,
    logger: logging.Logger,
    profile_path: str,
) -> None:
    """Ensure SQLite database exists, create if needed.

    Args:
        database_path: Path to SQLite database
        logger: Logger object
        profile_path: Path to rMLST TSV profile file
    """
    logger.info('Creating SQLite database: %s', database_path)

    # Create database connection
    conn = sqlite3.connect(database_path)

    # Read and convert profiles
    logger.debug('Reading profiles from TSV')
    profiles = pd.read_csv(
        profile_path,
        sep='\t'
    )

    # Create table
    logger.debug('Creating profiles table')
    profiles.to_sql(
        'profiles',
        conn,
        if_exists='replace',
        index=False
    )

    # Create indexes
    logger.debug('Creating database indexes')
    conn.execute('CREATE INDEX idx_rst ON profiles(rST)')
    conn.execute('CREATE INDEX idx_species ON profiles(species)')

    conn.close()
    logger.info('Database creation complete')


def _get_rmlst_profiles(
    *,
    conn: sqlite3.Connection,
    logger: logging.Logger,
    sample: CustomBox
) -> CustomBox:
    """Query rMLST profiles matching sample alleles.

    Args:
        conn: SQLite database connection
        logger: Logger object
        sample: Sample metadata object

    Returns:
        Updated sample object with matching profiles
    """
    # Build query conditions for each locus
    conditions = []
    for locus, alleles in sample.rmlst.allele_calls.items():
        allele_list = ','.join(alleles)
        conditions.append(f"{locus} IN ({allele_list})")

    # Create and execute query
    query = f"""
    SELECT rST, genus, species, subspecies
    FROM profiles
    WHERE {' AND '.join(conditions)}
    """

    # Fetch results
    matches = conn.execute(query).fetchall()

    # Store results in sample metadata
    sample.rmlst.matches = {
        'rSTs': [m[0] for m in matches],
        'genus': list(set(m[1] for m in matches)),
        'species': list(set(m[2] for m in matches))
    }

    logger.debug(
        'Sample %s matches rSTs: %s',
        sample.name,
        sample.rmlst.matches['rSTs']
    )

    return sample


def _set_metadata_attributes(
    *,  # Enforce keyword arguments
    reference_file_path: str,
    sample: CustomBox
) -> List[CustomBox]:
    """
    Create attributes for the metadata objects

    Args:
        reference_file_path (str): Path to the reference database
        sample (CustomBox): Metadata object for the sample

    Returns:
        CustomBox: Updated metadata object
    """
    # Set the rMLST database
    rmlst_database = os.path.join(
        reference_file_path,
        'rmlst',
        'rmlst'
    )

    # Create t~he rMLST attribute
    if not sample.key_exists('rmlst'):
        sample.rmlst = CustomBox()

    # Create the rMLST directory attribute
    sample.rmlst.output_dir = os.path.join(
        sample.general.output_directory,
        'rmlst'
    )

    # Create the output directory
    os.makedirs(sample.rmlst.output_dir, exist_ok=True)

    # Set the name of the KMA rMLST report
    sample.rmlst.report = os.path.join(
        sample.rmlst.output_dir,
        'rmlst.res'
    )

    # Set the name of the outputs
    sample.rmlst.outputs = os.path.join(
        sample.rmlst.output_dir,
        'rmlst'
    )
    # Set the KMA rMLST command
    sample.commands.rmlst = (
        f'kma -ipe {' '.join(sample.general.trimmed_corrected_fastq_files)} '
        f'-o {sample.rmlst.outputs} '
        f'-t_db {rmlst_database} '
        f'-ID 100'
    )

    return sample


def _parse_kma_results(
    *,  # Enforce keyword arguments
    cutoff: float = 100.0,
    logger: logging.Logger,
    min_depth: float = 5.0,
    sample: CustomBox,
) -> CustomBox:
    """Parse KMA rMLST results, grouping multiple alleles by locus.

    Args:
        cutoff: Minimum percent identity threshold
        logger: Logger object
        min_depth: Minimum depth of coverage threshold
        sample: CustomBox object containing sample metadata

    Returns:
        Updated sample object with parsed results grouped by locus
    """
    # Initialize results structure
    sample.rmlst.results = {}

    with open(sample.rmlst.report, 'r', encoding='utf-8') as res_file:
        next(res_file)
        for line in res_file:
            fields = line.strip().split('\t')

            template = fields[0]
            locus = template.split('_')[0]
            allele = template.split('_')[1]
            identity = float(fields[4])
            depth = float(fields[8])

            if identity >= cutoff and depth >= min_depth:
                # Initialize locus dict if needed
                if locus not in sample.rmlst.results:
                    sample.rmlst.results[locus] = {}

                # Store result data
                sample.rmlst.results[locus][allele] = {
                    'identity': identity,
                    'depth': depth,
                    'template_coverage': float(fields[5]),
                    'query_coverage': float(fields[7])
                }

    # Get best allele call for each locus
    sample.rmlst.allele_calls = _get_allele_calls(sample=sample)

    # Log the results
    logger.debug(
        'Sample: %s, results: %s',
        sample.name, sample.rmlst.allele_calls
    )

    return sample


def _get_allele_calls(
    # *, # Enforce keyword arguments
    sample: CustomBox,
) -> dict[str, list[str]]:
    """Extract allele calls meeting depth criteria for each locus.

    Args:
        sample: CustomBox object with parsed rMLST results

    Returns:
        Dictionary mapping locus IDs to lists of matching allele numbers
    """
    allele_calls = {}

    for locus, alleles in sample.rmlst.results.items():

        # Store the alleles in the allele_calls dictionary
        allele_calls[locus] = list(alleles.keys())

    return allele_calls
