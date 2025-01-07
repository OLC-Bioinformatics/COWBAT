#!/usr/env/python3

"""
Perform rMLST analyses on raw reads using KMA
"""

# Standard imports
from concurrent.futures import as_completed, ThreadPoolExecutor
import logging
import os
import sqlite3
from typing import (
    Dict,
    List,
    Tuple
)

# Third-party imports
from Bio import SeqIO
from olctools.accessoryFunctions.accessoryFunctions import run_subprocess
from olctools.accessoryFunctions.metadata import CustomBox

# Local imports
from cowbat.methods import write_to_log_files
from cowbat.rmlst_sqlite import (
    _get_rmlst_profiles,
    _sqlite_database
)

__author__ = 'adamkoziol'


def rmlst(
    *,  # Enforce keyword arguments
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    reference_file_path: str,
    report_path: str,
    threads: int
) -> List[CustomBox]:
    """
    Perform rMLST analyses on raw reads using KMA

    Args:
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger for recording information.
        metadata (List[CustomBox]): List of metadata objects for the samples.
        reference_file_path (str): Path to the reference database.
        report_path (str): Path to save the report.
        threads (int): Number of threads to use.

    Returns:
        List[CustomBox]: Updated metadata after all processing steps.
    """
    # Ensure that the SQLite reference database exists
    rmlst_database, profile = _sqlite_database(
        logger=logger,
        reference_file_path=reference_file_path
    )

    # Get rMLST gene names
    rmlst_genes = _get_rmlst_genes(profile_path=profile)

    # Set the path to the rMLST database
    rmlst_database_path = os.path.join(
        reference_file_path,
        'rmlst',
    )

    # Run the KMA analyses on the samples
    for sample in metadata:
        # Create the necessary attributes for the metadata objects
        sample = _set_metadata_attributes(
            reference_file_path=reference_file_path,
            sample=sample,
            threads=threads
        )

        # Run the KMA command on the filtered database
        if not os.path.isfile(sample.rmlst.report):
            # Log the command
            logger.debug(sample.commands.rmlst)
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
        # Ensure that all the rMLST genes are present in the sample
        unmapped, imperfect_matches = _identify_unmapped_genes(
            logger=logger,
            rmlst_genes=rmlst_genes,
            sample=sample
        )

        # If there are unmapped genes, filter the database
        if unmapped:
            best_alleles = _gene_specific_rmlst_database(
                imperfect_matches=imperfect_matches,
                log_file=log_file,
                logger=logger,
                sample=sample,
                rmlst_database_path=rmlst_database_path,
                threads=threads,
                unmapped_genes=unmapped
            )

            # Update the report with the best allele calls
            _update_report(
                best_alleles=best_alleles,
                logger=logger,
                sample=sample
            )

    # Load the SQLite database
    conn = sqlite3.connect(rmlst_database)

    # Set the name of the report
    report = os.path.join(
        report_path,
        'rmlst.tsv'
    )

    # Delete the report if it exists
    if os.path.isfile(report):
        os.remove(report)

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
            report=report,
            sample=sample
        )

    # Close database connection
    conn.close()

    return metadata


def _get_rmlst_genes(
    *,  # Enforce keyword arguments
    profile_path: str
) -> set[str]:
    """
    Extract rMLST gene names from profile file header

    Args:
        profile_path: Path to profile.txt

    Returns:
        Set of BACT* gene names
    """
    with open(profile_path, 'r', encoding='utf-8') as rmlst_database:
        header = rmlst_database.readline().strip().split('\t')
    return {col for col in header if col.startswith('BACT')}


def _set_metadata_attributes(
    *,  # Enforce keyword arguments
    reference_file_path: str,
    sample: CustomBox,
    threads: int
) -> List[CustomBox]:
    """
    Create attributes for the metadata objects

    Args:
        reference_file_path (str): Path to the reference database
        sample (CustomBox): Metadata object for the sample
        threads (int): Number of threads to use

    Returns:
        CustomBox: Updated metadata object
    """
    # Set the rMLST database
    rmlst_database = os.path.join(
        reference_file_path,
        'rmlst',
        'rmlst'
    )

    # Create the rMLST attribute
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

    # Set the name of the KMA rMLST analyses
    sample.rmlst.outputs = os.path.join(
        sample.rmlst.output_dir,
        'rmlst'
    )

    # Set the KMA rMLST system call
    sample.commands.rmlst = (
        f'kma -ipe {" ".join(sample.general.trimmed_corrected_fastq_files)} '
        f'-o {sample.rmlst.outputs} '
        f'-t_db {rmlst_database} '
        f'-t {threads} '    # Threading
        '-ID 50 '           # Low threshold for imperfect matches
        '-mem_mode '        # Memory mode
        '-a '               # Output all matches
    )
    return sample


def _identify_unmapped_genes(
    *,  # Enforce keyword arguments
    logger: logging.Logger,
    rmlst_genes: set[str],
    sample: CustomBox
) -> Tuple[set[str], dict[str, list[str]]]:
    """
    Identify genes that are not present in the rMLST database

    Args:
        logger (logging.Logger): Logger object
        rmlst_genes (set[str]): Set of rMLST gene names
        sample (CustomBox): Metadata object for the sample

    Returns:
        CustomBox: Updated metadata object
    """
    # Create a list of all templates in the counts file
    genes = set()

    # Create a dictionary to store all imperfect matches
    imperfect_matches = {}

    # Parse the KMA results
    with open(sample.rmlst.report, 'r', encoding='utf-8') as res_file:
        # Skip header
        next(res_file)
        # Iterate over the lines in the file
        for line in res_file:
            fields = line.strip().split('\t')

            # Extract the allele
            allele = fields[0]

            # Extract the gene name
            gene = allele.split('_')[0]

            # Extract the identity
            identity = float(fields[4])

            # Filter the imperfect matches
            if identity < 100.0:
                if gene not in imperfect_matches:
                    imperfect_matches[gene] = []
                imperfect_matches[gene].append(allele)
            else:
                # Add the gene to the set
                genes.add(gene)

    # Use a set comparison to identify genes that do not have a match in the
    # KMA outputs
    unmapped_genes = rmlst_genes - genes

    # Log the results
    logger.debug(
        'Sample: %s, unmapped genes: %s',
        sample.name, unmapped_genes
    )

    return unmapped_genes, imperfect_matches


def _gene_specific_rmlst_database(
    *,  # Enforce keyword arguments
    imperfect_matches: dict[str, list[str]],
    log_file: str,
    logger: logging.Logger,
    rmlst_database_path: str,
    sample: CustomBox,
    threads: int,
    unmapped_genes: set[str],
) -> Dict[str, str]:
    """
    Create a FASTA file of the unmapped

    Args:
        imperfect_matches (dict[str, list[str]]): Dictionary of gene names
            and lists of imperfect
        log_file (str): Path to the
        logger (logging.Logger): Logger object
        rmlst_database_path (str): Path to the rMLST database
        sample (CustomBox): Metadata object for the sample
        threads (int): Number of threads
        unmapped_genes (set[str]): Set of unmapped gene names
        num_alleles (int): Maximum number of alleles to include in
            sub-databases

    Returns:
        Dict[str, str]: Dictionary of gene names and best allele calls
    """

    best_alleles_dict = {}
    # Iterate over the unmapped genes
    for gene in unmapped_genes:
        # Create the FASTA file
        gene_database_path, gene_alleles = _link_rmlst_gene_file(
            gene=gene,
            logger=logger,
            sample=sample,
            rmlst_database_path=rmlst_database_path
        )

        # Set the KMA database name
        gene_database_name = os.path.join(
            gene_database_path,
            gene
        )

        # Create a KMA index command for the filtered database
        kma_gene_index_command = (
            f'kma index -i {gene_alleles} -o {gene_database_name}'
        )

        # Log the command
        logger.debug(kma_gene_index_command)

        # Run the KMA index command
        if not os.path.isfile(f'{gene_database_name}.name'):
            out, err = run_subprocess(
                command=kma_gene_index_command
            )

            # Write the command and the outputs to the log files
            write_to_log_files(
                command=kma_gene_index_command,
                err=err,
                log_file=log_file,
                logger=logger,
                out=out,
                program='rmlst',
                sample=sample
            )

        # Create the command to run KMA on the gene database
        kma_gene_command = (
            'kma '
            f'-ipe {" ".join(sample.general.trimmed_corrected_fastq_files)} '
            f'-o {gene_database_name} '
            f'-t_db {gene_database_name} '
            f'-t {threads} '    # Threading
            '-mem_mode '        # Memory mode
            '-sasm '            # Skip alignment
            '-ck'               # Count k-mers over pseudo alignment
        )

        # Log the command
        logger.debug(kma_gene_command)

        # Set the name of the report
        gene_specific_report = os.path.join(
            gene_database_path,
            f'{gene}.res'
        )

        # Run the KMA command on the gene database
        if not os.path.isfile(gene_specific_report):
            out, err = run_subprocess(
                command=kma_gene_command
            )

            # Write the command and the outputs to the log files
            write_to_log_files(
                command=kma_gene_command,
                err=err,
                log_file=log_file,
                logger=logger,
                out=out,
                program='rmlst',
                sample=sample
            )

        # Parse the KMA results to extract all matches exceeding the p-value
        # threshold
        best_alleles, best_line = _filter_database(
            gene=gene,
            gene_alleles=gene_alleles,
            gene_database_path=gene_database_path,
            gene_specific_report=gene_specific_report,
            imperfect_matches=imperfect_matches,
            log_file=log_file,
            logger=logger,
            sample=sample,
            threads=threads
        )

        best_alleles_dict[best_alleles] = best_line
    return best_alleles_dict


def _link_rmlst_gene_file(
    *,  # Enforce keyword arguments
    gene: str,
    logger: logging.Logger,
    rmlst_database_path: str,
    sample: CustomBox
) -> Tuple[str, str]:
    """
    Create a symlink to the rMLST gene file in the sample directory

    Args:
        gene (str): Name of the gene
        logger (logging.Logger): Logger object
        rmlst_database_path (str): Path to the rMLST database
        sample (CustomBox): Metadata object for the sample

    Returns:
        Tuple[str, str]: Path to the gene database and the filtered FASTA file
    """

    # Set the name and path of the gene file
    gene_file = os.path.join(
        rmlst_database_path,
        f'{gene}.tfa'
    )

    # Set the path to the gene database
    gene_database_path = os.path.join(
        sample.rmlst.output_dir,
        gene
    )

    # Create the directory
    os.makedirs(gene_database_path, exist_ok=True)

    # Set the path to the filtered rMLST database
    gene_alleles = os.path.join(
        gene_database_path,
        f'{gene}.fasta'
    )

    # Create a symlink to the gene file
    if not os.path.isfile(gene_alleles):
        os.symlink(
            gene_file,
            gene_alleles
        )
        # Log the symlink creation
        logger.debug('Created symlink to %s', gene_alleles)

    return gene_database_path, gene_alleles


def _filter_database(
    *,  # Enforce keyword arguments
    gene: str,
    gene_alleles: str,
    gene_database_path: str,
    gene_specific_report: str,
    imperfect_matches: dict[str, list],
    log_file: str,
    logger: logging.Logger,
    sample: CustomBox,
    threads: int,
    max_alleles: int = 1,
    p_value_threshold: float = 1.0e-20
) -> Tuple[str, str]:
    """
    Filter the rMLST database to only include loci present in the sample

    Args:
        gene (str): Name of the unmapped gene
        gene_alleles (str): Path to the gene-specific FASTA file
        gene_database_path (str): Path to the gene database
        gene_specific_report (str): Path to the gene-specific KMA report
        imperfect_matches (dict[str, list[str]]): Dictionary of gene names
            and lists of imperfect matches
        log_file (str): Path to the log file
        logger (logging.Logger): Logger object
        sample (CustomBox): Metadata object for the sample
        threads (int): Number of threads to use
        max_alleles (int): Maximum number of alleles to include in the
            sub-databases. Default is 1
        p_value_threshold (float): P-value threshold for filtering alleles

    Returns:
        Tuple[str, str]: Path to the gene database and the filtered FASTA file
    """

    # Parse the KMA report and extract all the alleles that exceed the p-value
    # threshold
    alleles = set()

    with open(gene_specific_report, 'r', encoding='utf-8') as res_file:
        # Skip header
        next(res_file)
        # Iterate over the lines in the file
        for line in res_file:
            fields = line.strip().split('\t')
            template = fields[0]
            p_value = float(fields[-1])
            if p_value < p_value_threshold:
                # Ensure that the allele is not an imperfect match
                try:
                    if template not in imperfect_matches[gene]:
                        alleles.add(template)
                    else:
                        logger.debug('Imperfect match: %s', template)
                except KeyError:
                    alleles.add(template)
    sequences = []

    # Set the path to the filtered rMLST database
    gene_database_path = gene_database_path + '_filtered'
    os.makedirs(gene_database_path, exist_ok=True)

    # Extract the sequence of the alleles from the gene-specific FASTA file
    with open(gene_alleles, 'r', encoding='utf-8') as fasta_file:
        records = SeqIO.parse(fasta_file, 'fasta')
        for record in records:
            if record.id in alleles:
                sequences.append(record)

    # Determine the number of sub-databases to create
    num_sub_databases = len(alleles) // max_alleles + 1

    # Initialize lists to store the KMA commands and report files
    kma_commands = []
    report_files = []

    # Iterate over the sub-databases
    for database_number in range(num_sub_databases):
        # Set the path to the sub-database
        sub_database = os.path.join(
            gene_database_path,
            f'{gene}_{database_number}.fasta'
        )

        # Log the sub-database creation
        logger.debug('Creating sub-database: %s', sub_database)

        # Write the alleles to the sub-database file
        with open(sub_database, 'w', encoding='utf-8') as sub_file:
            SeqIO.write(
                sequences[
                    database_number * max_alleles:(
                        database_number + 1
                    ) * max_alleles
                ], sub_file, 'fasta'
            )

        # Index the sub-database
        sub_database_index = os.path.join(
            gene_database_path,
            f'{gene}_{database_number}'
        )

        # Create the KMA index command for the sub-database
        kma_sub_index_command = (
            f'kma index -i {sub_database} -o {sub_database_index}'
        )

        # Log the command
        logger.debug(kma_sub_index_command)

        # Set the name of the sub-database index output file
        sub_database_index_output = os.path.join(
            gene_database_path,
            f'{sub_database_index}.name'
        )

        if not os.path.isfile(sub_database_index_output):
            # Run the KMA index command
            out, err = run_subprocess(
                command=kma_sub_index_command
            )

            # Write the command and the outputs to the log files
            write_to_log_files(
                command=kma_sub_index_command,
                err=err,
                log_file=log_file,
                logger=logger,
                out=out,
                program='rmlst',
                sample=sample
            )

        # Set the name of the sub-database report
        sub_database_report = os.path.join(
            gene_database_path,
            f'{sub_database_index}.res'
        )

        # Create the KMA command for the sub-database analysis
        kma_sub_command = (
            'kma '
            f'-ipe {" ".join(sample.general.trimmed_corrected_fastq_files)} '
            f'-o {sub_database_index} '
            f'-t_db {sub_database_index} '
            # f'-t {threads} '   # Threading
            '-ID 99 '          # Allow imperfect matches
            '-mem_mode '       # Memory mode
        )

        # Log the command
        logger.debug(kma_sub_command)

        # Store commands for parallel execution
        kma_commands.append(kma_sub_command)
        report_files.append(sub_database_report)

    # Parse the KMA results to extract the best allele call
    best_allele = None
    best_allele_identity = 0.0
    best_line = None

    # Run KMA commands in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                _run_kma_search,
                kma_command=cmd,
            )
            for cmd, report in zip(kma_commands, report_files)
            if not os.path.isfile(report)
        ]
        for future in as_completed(futures):
            future.result()

    # Process each sub-database
    for database_number in range(num_sub_databases):

        # Index the sub-database
        sub_database_index = os.path.join(
            gene_database_path,
            f'{gene}_{database_number}'
        )
        # Set the name of the sub-database report
        sub_database_report = os.path.join(
            gene_database_path,
            f'{sub_database_index}.res'
        )
        with open(sub_database_report, 'r', encoding='utf-8') as res_file:
            try:
                next(res_file)
                for line in res_file:
                    fields = line.strip().split('\t')
                    template = fields[0]
                    identity = float(fields[4])
                    if identity > best_allele_identity:
                        best_allele = template
                        best_allele_identity = identity
                        best_line = line
                        logger.debug(
                            'Allele: %s, identity: %.1f%%', template, identity
                        )
            except StopIteration:
                pass

    # Log the best allele call
    logger.debug(
        'Best allele call: %s, identity: %.1f%%',
        best_allele, best_allele_identity
    )

    return best_allele, best_line


def _run_kma_search(
    *,
    kma_command: str,
) -> None:
    """
    Run KMA search and handle logging

    Args:
        kma_command (str): KMA search command
    """
    run_subprocess(command=kma_command)


def _update_report(
    *,  # Enforce keyword arguments
    best_alleles: dict[str, str],
    logger: logging.Logger,
    sample: CustomBox
):
    """
    Update the KMA report with the best allele calls

    Args:
        best_alleles (dict[str, str]): Dictionary of gene names and
            best allele calls
        logger (logging.Logger): Logger object
        sample (CustomBox): Metadata object for the sample
    """
    # Read existing lines
    with open(sample.rmlst.report, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Parse lines and find insertion points
    updated_lines = [lines[0]]  # Keep header
    current_gene = None

    for line in lines[1:]:  # Skip header
        fields = line.strip().split('\t')
        gene = fields[0].split('_')[0]  # Extract BACT number

        # Check if we need to insert any new results before this gene
        if current_gene != gene:
            for allele, result in best_alleles.items():
                if (
                    allele.split('_')[0] < gene
                    and allele not in [
                        line.split('\t')[0] for line in updated_lines
                    ]
                ):
                    logger.debug('Inserting result: %s', result)
                    updated_lines.append(result)

        # Add current line
        updated_lines.append(line)
        current_gene = gene

    # Add any remaining results at end
    for allele, result in best_alleles.items():
        if allele not in [line.split('\t')[0] for line in updated_lines]:
            updated_lines.append(result + '\n')

    # Write back to file
    with open(sample.rmlst.report, 'w', encoding='utf-8') as updated_file:
        updated_file.writelines(updated_lines)


def _parse_kma_results(
    *,  # Enforce keyword arguments
    cutoff: float = 100.0,
    logger: logging.Logger,
    sample: CustomBox,
) -> CustomBox:
    """Parse KMA rMLST results, grouping multiple alleles by locus.

    Args:
        cutoff: Minimum percent identity threshold
        logger: Logger object
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

            if identity >= cutoff:
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

    # Iterate over loci and alleles
    for locus, alleles in sample.rmlst.results.items():

        # Store the alleles in the allele_calls dictionary
        allele_calls[locus] = list(alleles.keys())

    return allele_calls
