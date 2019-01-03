#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import SetupLogging
from validator_helper import validate
from argparse import ArgumentParser
import logging


def validate_cowbat(reference_report, test_report, assembly_typer=False):
    column_list = list()
    column_list.append(validate.Column(name='Genus'))
    column_list.append(validate.Column(name='N50', column_type='Range', acceptable_range=40000))
    column_list.append(validate.Column(name='NumContigs', column_type='Range', acceptable_range=50))
    column_list.append(validate.Column(name='TotalLength', column_type='Range', acceptable_range=100000))
    column_list.append(validate.Column(name='PercentGC', column_type='Range', acceptable_range=2))
    column_list.append(validate.Column(name='rMLST_Result'))
    column_list.append(validate.Column(name='MLST_Result'))
    column_list.append(validate.Column(name='E_coli_Serotype'))
    column_list.append(validate.Column(name='SISTR_serovar_cgMLST'))
    column_list.append(validate.Column(name='GeneSeekr_Profile'))
    column_list.append(validate.Column(name='Vtyper_Profile'))
    column_list.append(validate.Column(name='AMR_Profile'))
    column_list.append(validate.Column(name='PlasmidProfile'))
    column_list.append(validate.Column(name='TotalPredictedGenes', column_type='Range', acceptable_range=50))
    if not assembly_typer:
        column_list.append(validate.Column(name='SamplePurity'))
        column_list.append(validate.Column(name='MeanInsertSize', column_type='Range', acceptable_range=50))
        column_list.append(validate.Column(name='InsertSizeSTD', column_type='Range', acceptable_range=50))
        column_list.append(validate.Column(name='AverageCoverageDepth', column_type='Range', acceptable_range=10))
        column_list.append(validate.Column(name='CoverageDepthSTD', column_type='Range', acceptable_range=10))

    validator = validate.Validator(reference_csv=reference_report,
                                   test_csv=test_report,
                                   column_list=column_list,
                                   identifying_column='SeqID')

    validation_list = [validator.same_columns_in_ref_and_test(),
                       validator.all_test_columns_in_ref_and_test(),
                       validator.check_samples_present(),
                       validator.check_columns_match()]
    if False in validation_list:
        return False
    else:
        return True


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Run integration tests on COWBAT pipeline')
    parser.add_argument('-r', '--reference_csv',
                        required=True,
                        help='Path to reference CSV with acceptable data.')
    parser.add_argument('-t', '--test_csv',
                        required=True,
                        help='Path to test CSV with data to be evaluated.')
    parser.add_argument('-a', '--assembly',
                        action='store_true',
                        default=False,
                        help='The assembly typing pipeline was used to process the run, rather than full COWBAT')
    # Get the arguments into an object
    args = parser.parse_args()
    # Pretty logging!
    SetupLogging()
    # Test the reports.
    successful_validate = validate_cowbat(reference_report=args.reference_csv,
                                          test_report=args.test_csv,
                                          assembly_typer=args.assembly)
    if successful_validate:
        logging.info('COWBAT successfully validated! :D')
    else:
        logging.error('COWBAT not successfully validated.')
