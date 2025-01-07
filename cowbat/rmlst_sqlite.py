#!/usr/bin/env python3

"""
Methods for rMLST typing using SQLite database
"""

# Standard imports
import logging
import os
import sqlite3
from typing import Tuple

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import CustomBox
import pandas as pd

__author__ = 'adamkoziol'


def _sqlite_database(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    reference_file_path: str
) -> Tuple[str, str]:
    """
    Ensure that the SQLite reference database exists

    Args:
        logger (logging.Logger): Logger object
        reference_file_path (str): Path to the reference database

    Returns:
        Tuple[str, str]: Path to the SQLite database and the rMLST profile file
    """
    # Set the path to the SQLite database
    sqlite_database = os.path.join(
        reference_file_path,
        'rmlst',
        'rmlst.sqlite'
    )

    # Set the path to the rMLST profile file
    profile_path = os.path.join(
            reference_file_path,
            'rmlst',
            'profile.txt'
        )

    # Ensure that the SQLite database exists
    if not os.path.isfile(sqlite_database):

        # Create the SQLite database
        _create_rmlst_database(
            database_path=sqlite_database,
            logger=logger,
            profile_path=profile_path
        )
    return sqlite_database, profile_path


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
    *,  # Enforce keyword arguments
    conn: sqlite3.Connection,
    logger: logging.Logger,
    report: str,
    sample: CustomBox
) -> CustomBox:
    """
    Query rMLST profiles matching sample alleles.

    Args:
        conn: SQLite database connection
        logger: Logger object
        report: Path to the KMA rMLST report
        sample: Sample metadata object

    Returns:
        Updated sample object with matching profiles
    """
    cursor = conn.cursor()

    # Get column names
    cursor.execute("PRAGMA table_info(profiles)")
    db_columns = [col[1] for col in cursor.fetchall()]

    # Define output column order
    output_columns = [
        'Strain', 'Genus', 'rMLST_genus', 'species', 'subspecies', 'lineage',
        'sublineage', 'other_designation', 'notes', 'SequenceType', 'Matches'
    ]
    # Add BACT genes to columns
    bact_columns = [col for col in db_columns if col.startswith('BACT')]
    output_columns.extend(bact_columns)

    # Create/check report file
    file_exists = os.path.isfile(report)
    mode = 'a' if file_exists else 'w'

    with open(report, mode, encoding='utf-8') as report:
        # Write header if new file
        if not file_exists:
            report.write('\t'.join(output_columns) + '\n')

        # Phase 1: Try perfect matches first
        conditions = []
        for locus, alleles in sample.rmlst.allele_calls.items():
            if len(alleles) == 1:
                conditions.append(
                    f"({locus} = '{alleles[0]}' OR {locus} = 'N')"
                )
            else:
                conditions.append(f"{locus} = 'N'")

        query = f"""
        SELECT *
        FROM profiles
        WHERE {' AND '.join(conditions)}
        """

        # Execute query
        cursor.execute(query)
        perfect_matches = cursor.fetchall()

        if perfect_matches:
            logger.debug("Found perfect matches!")
            matches = []
            for match in perfect_matches:
                profile_dict = dict(zip(db_columns, match))
                matches.append(
                    {
                        'rST': profile_dict['rST'],
                        'genus': profile_dict['genus'],
                        'species': profile_dict['species'],
                        'subspecies': profile_dict.get('subspecies', ''),
                        'lineage': profile_dict.get('lineage', ''),
                        'sublineage': profile_dict.get('sublineage', ''),
                        'other_designation': profile_dict.get(
                            'other_designation', ''
                        ),
                        'notes': profile_dict.get('notes', ''),
                        'matches': len(conditions),
                        'total': len(conditions),
                        'percent': 100.0,
                        'alleles': {
                            col: (val, sample.rmlst.allele_calls.get(col, []))
                            for col, val in zip(db_columns, match)
                        }
                    }
                )
        else:
            # Phase 2: Find best partial matches
            logger.debug("No perfect matches, finding best partial matches...")
            cursor.execute("SELECT * FROM profiles")
            profiles = cursor.fetchall()

            # Initialize best count and matches
            best_count = 0
            matches = []

            for profile in profiles:
                match_count = 0
                profile_dict = dict(zip(db_columns, profile))
                allele_mapping = {}

                # Count matches and track alleles
                for locus, alleles in sample.rmlst.allele_calls.items():
                    db_value = profile_dict[locus]
                    allele_mapping[locus] = (db_value, alleles)
                    if db_value == 'N' or any(
                        str(allele) == str(db_value) for allele in alleles
                    ):
                        match_count += 1

                # Only keep if it's the best match so far
                if match_count >= best_count:
                    if match_count > best_count:
                        matches = []  # Clear previous
                        best_count = match_count
                    matches.append({
                        'rST': profile_dict['rST'],
                        'genus': profile_dict['genus'],
                        'species': profile_dict['species'],
                        'subspecies': profile_dict.get('subspecies', ''),
                        'lineage': profile_dict.get('lineage', ''),
                        'sublineage': profile_dict.get('sublineage', ''),
                        'other_designation': profile_dict.get(
                            'other_designation', ''
                        ),
                        'notes': profile_dict.get('notes', ''),
                        'matches': match_count,
                        'total': len(sample.rmlst.allele_calls),
                        'percent': (
                            match_count / len(sample.rmlst.allele_calls)
                        ) * 100,
                        'alleles': allele_mapping
                    })

        # Write matching profiles
        if matches:
            for match in matches:
                # Build profile line in new column order
                logger.debug("Writing profile: %s", match)
                profile_line = [
                    str(sample.name),  # Strain
                    str(sample.general.reference_genus),  # Genus
                    str(match['genus']),  # rMLST_genus
                    str(match['species']),
                    str(match.get('subspecies', '')),  # Convert None to ''
                    str(match.get('lineage', '')),
                    str(match.get('sublineage', '')),
                    str(match.get('other_designation', '')),
                    str(match.get('notes', '')),
                    str(match['rST']),
                    str(match['matches'])
                ]

                # Add formatted alleles
                for col in bact_columns:
                    if col in match['alleles']:
                        expected, observed = match['alleles'][col]
                        profile_line.append(_format_allele_display(
                            observed=observed,
                            expected=expected
                        ))
                    else:
                        profile_line.append('N')

                report.write('\t'.join(profile_line) + '\n')

    # Store results
    sample.rmlst.matches = {
        'rSTs': [m['rST'] for m in matches],
        'genus': list(set(m['genus'] for m in matches)),
        'species': list(set(m['species'] for m in matches)),
        'match_details': matches
    }

    logger.debug(
        'Sample %s matches rSTs: %s (best match: %d/%d = %.1f%%)',
        sample.name,
        sample.rmlst.matches['rSTs'],
        matches[0]['matches'],
        matches[0]['total'],
        matches[0]['percent']
    )
    return sample


def _format_allele_display(
    *,  # Enforce keyword arguments
    observed: list,
    expected: str
) -> str:
    """
    Format allele display based on match type

    Args:
        observed (list): List of observed alleles
        expected (str): Expected allele

    Returns:
        str: Formatted allele display

    Cases:
    - Perfect match: "10"
    - Multiple alleles: "10 692 (N)"
    - Missing allele: "NA (45)"
    """
    if not observed:  # Missing allele
        return f"NA ({expected})"
    elif len(observed) > 1:  # Multiple alleles
        return f"{' '.join(map(str, observed))} (N)"
    elif str(observed[0]) != str(expected) and expected != 'N':  # Mismatch
        return f"{expected} ({observed[0]})"
    else:  # Perfect match or N
        return str(expected)
