#!/usr/bin/env python3

"""
Confindr contamination detection in raw reads
"""

# Standard imports
import logging
import os
from typing import (
    Any,
    Dict,
    List,
    Optional
)

# Third-party imports
import pandas as pd
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def contamination_finder(
    input_path: str,
    log_file: str,
    logger: logging.Logger,
    metadata: List[CustomBox],
    debug: bool = False,
    report_path: Optional[str] = None,
    threads: int = 12
) -> None:
    """
    Helper function to get ConFindr integrated into the assembly pipeline.

    Args:
        input_path (str): Path to the input directory.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.
        debug (bool): Whether to keep intermediate files for debugging.
        report_path (Optional[str]): Path to the report directory.
        threads (int): Number of threads to use.
    """
    logger.info('Calculating contamination in reads')

    # Set the path variables
    report_path = report_path if report_path else os.path.join(
        input_path, 'confindr')

    # Create the report path if required
    os.makedirs(report_path, exist_ok=True)

    confindr_report = os.path.join(
        input_path, 'confindr', 'confindr_report.csv')
    pipeline_report = os.path.join(report_path, 'confindr_report.tsv')

    # Only proceed if the ConFindr report doesn't exist
    if not os.path.isfile(confindr_report):
        run_confindr(
            debug=debug,
            log_file=log_file,
            logger=logger,
            input_path=input_path,
            reportpath=report_path,
            threads=threads,
        )

    # Load the ConFindr report into a dictionary using pandas
    confindr_results = load_confindr_report(confindr_report=confindr_report)

    # Process the ConFindr results for each sample
    metadata = process_confindr_results(
        confindr_results=confindr_results,
        logger=logger,
        metadata=metadata,
    )

    # Write the processed results to the pipeline report
    write_pipeline_report(
        metadata=metadata,
        pipeline_report=pipeline_report
    )

    return metadata


def run_confindr(
    debug: bool,
    log_file: str,
    input_path: str,
    logger: logging.Logger,
    reportpath: str,
    threads: int
) -> None:
    """
    Run ConFindr to detect contamination in reads.

    Args:
        debug (bool): Whether to keep intermediate files for debugging.
        log_file (str): Path to the log file.
        logger (logging.Logger): Logger object.
        input_path (str): Path to the input directory.
        reportpath (str): Path to the report directory.
        threads (int): Number of threads to use.
    """
    os.makedirs(reportpath, exist_ok=True)
    logger.debug('Running ConFindr for input directory: %s', input_path)

    # Create the necessary paths
    output_path = os.path.join(input_path, "confindr")

    # Construct the ConFindr system call
    system_call = (
        f'confindr.py -i {input_path} -o {output_path} -t {threads} --rmlst'
    )

    # Keep the files if the debug option is specified
    if debug:
        system_call += ' -k'

    # Run the ConFindr system call
    out, err = run_subprocess(system_call)
    write_to_log_file(
        out=system_call,
        err=system_call,
        log_file=log_file
    )
    write_to_log_file(
        out=out,
        err=err,
        log_file=log_file,
        sample_log=None,
        sample_err=None
    )
    logger.info('Contamination detection complete!')


def load_confindr_report(confindr_report: str) -> Dict[str, Dict[str, Any]]:
    """
    Load the ConFindr report into a dictionary using pandas.

    Args:
        confindr_report (str): Path to the ConFindr report.

    Returns:
        Dict[str, Dict[str, Any]]: Dictionary containing ConFindr results.
    """
    return pd.read_csv(confindr_report, index_col=0).T.to_dict()


def process_confindr_results(
    confindr_results: Dict[str, Dict[str, Any]],
    logger: logging.Logger,
    metadata: List[CustomBox],
) -> None:
    """
    Process the ConFindr results for each sample.

    Args:
        confindr_results (Dict[str, Dict[str, Any]]): Dictionary containing
        ConFindr results.
        logging (logging.Logger): Logger object.
        metadata (List[CustomBox]): List of metadata sample objects.
    """
    if not metadata or not confindr_results:
        logger.warning('No metadata or ConFindr results to process.')
        return

    for sample in metadata:
        logger.debug(
            'Processing ConFindr results for sample: %s',
            sample.name)
        sample.confindr = CustomBox()

        # Iterate through the dictionary to find the outputs for each sample
        for line in confindr_results:
            if sample.name not in line:
                continue

            # Set the values using the appropriate keys as the attributes
            sample.confindr.genus = confindr_results[line]['Genus'] \
                if isinstance(confindr_results[line]['Genus'], str) \
                else 'ND'
            sample.confindr.num_contaminated_snvs = \
                confindr_results[line]['NumContamSNVs']
            sample.confindr.contamination_status = confindr_results[line][
                'ContamStatus']

            # Handle percent contamination calculations
            sample.confindr.percent_contamination = get_confindr_value(
                confindr_results[line], 'PercentContam', 0
            )
            sample.confindr.percent_contamination_std = get_confindr_value(
                confindr_results[line], 'PercentContamStandardDeviation', 0
            )

            if sample.confindr.contamination_status == 'True':
                sample.confindr.contamination_status = 'Contaminated'
            elif sample.confindr.contamination_status == 'False':
                sample.confindr.contamination_status = 'Clean'

    return metadata


def get_confindr_value(
        result_line: Dict[str, Any],
        key: str, default: Any) -> Any:
    """
    Get a value from the ConFindr result line, handling missing keys.

    Args:
        result_line (Dict[str, Any]): Dictionary containing a line of ConFindr
        results.
        key (str): Key to look up in the result line.
        default (Any): Default value to return if the key is missing or invalid

    Returns:
        Any: Value from the result line or the default value.
    """
    try:
        value = result_line[key]
        return value if str(value) != 'nan' else default
    except KeyError:
        return default


def write_pipeline_report(
        metadata: List[CustomBox],
        pipeline_report: str) -> None:
    """
    Write the processed ConFindr results to the pipeline report.

    Args:
        metadata (List[CustomBox]): List of metadata sample objects.
        pipeline_report (str): Path to the pipeline report.
    """
    with open(pipeline_report, 'w', encoding='utf-8') as tsv:
        data = 'Strain,Genus,NumContamSNVs,ContamStatus,PercentContam,' \
            'PercentContamSTD\n'
        for sample in metadata:
            data += (
                f'{sample.name}\t'
                f'{sample.confindr.genus}\t'
                f'{sample.confindr.num_contaminated_snvs}\t'
                f'{sample.confindr.contamination_status}\t'
                f'{sample.confindr.percent_contamination}\t'
                f'{sample.confindr.percent_contamination_std}\n'
            )
        tsv.write(data)
