#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import SetupLogging
from validator_helper import validate
from argparse import ArgumentParser
import logging
import os


class ValidateCowbat(object):

    def validate_cowbat(self):
        # combinedMetadata
        logging.info('Validating combinedMetadata.csv')
        self.validate_combined_metadata(reference_report=os.path.join(self.reference_folder, 'combinedMetadata.csv'),
                                        test_report=os.path.join(self.test_folder, 'combinedMetadata.csv'),
                                        assembly_typer=self.assembly_typer)
        logging.info('Double validating combinedMetadata.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'combinedMetadata.csv'),
                             test_report=os.path.join(self.test_folder, 'combinedMetadata.csv'),
                             columns_to_exclude=['SeqID', 'AssemblyDate'],
                             identifying_column='SeqID')
        # AMRSummary
        logging.info('Validating amr_summary.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'amr_summary.csv'),
                             test_report=os.path.join(self.test_folder, 'amr_summary.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # ConFindr
        logging.info('Validating confindr_report.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'confindr_report.csv'),
                             test_report=os.path.join(self.test_folder, 'confindr_report.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # CoreGenome
        logging.info('Validating coregenome.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'coregenome.csv'),
                             test_report=os.path.join(self.test_folder, 'coregenome.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # Escherichia core genome
        logging.info('Validating Escherichia_core.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'Escherichia_core.csv'),
                             test_report=os.path.join(self.test_folder, 'Escherichia_core.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # GDCS
        logging.info('Validating GDCS.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'GDCS.csv'),
                             test_report=os.path.join(self.test_folder, 'GDCS.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # GeneSippr
        logging.info('Validating genesippr.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'genesippr.csv'),
                             test_report=os.path.join(self.test_folder, 'genesippr.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # Legacy V-typer
        logging.info('Validating legacy_vtyper.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'legacy_vtyper.csv'),
                             test_report=os.path.join(self.test_folder, 'legacy_vtyper.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # MASH
        logging.info('Validating mash.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mash.csv'),
                             test_report=os.path.join(self.test_folder, 'mash.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # MLST Bacillus
        logging.info('Validating mlst_Bacillus.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Bacillus.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Bacillus.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MLST Campylobacter
        logging.info('Validating mlst_Campylobacter.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Campylobacter.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Campylobacter.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MLST Escherichia
        logging.info('Validating mlst_Escherichia.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Escherichia.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Escherichia.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MLST Listeria
        logging.info('Validating mlst_Listeria.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Listeria.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Listeria.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MLST Salmonella
        logging.info('Validating mlst_Salmonella.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Salmonella.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Salmonella.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MLST Vibrio
        logging.info('Validating mlst_Vibrio.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mlst_Vibrio.csv'),
                             test_report=os.path.join(self.test_folder, 'mlst_Vibrio.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # MOB-Recon
        logging.info('Validating mob_recon_summary.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'mob_recon_summary.csv'),
                             test_report=os.path.join(self.test_folder, 'mob_recon_summary.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # Prophages
        logging.info('Validating prophages.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'prophages.csv'),
                             test_report=os.path.join(self.test_folder, 'prophages.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # QAML report
        logging.info('Validating QAMLReport.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'QAMLReport.csv'),
                             test_report=os.path.join(self.test_folder, 'QAMLReport.csv'),
                             columns_to_exclude=['Sample'],
                             identifying_column='Sample')
        # ResFinder
        logging.info('Validating resfinder.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'resfinder.csv'),
                             test_report=os.path.join(self.test_folder, 'resfinder.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # rMLST
        logging.info('Validating rmlst.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'rmlst.csv'),
                             test_report=os.path.join(self.test_folder, 'rmlst.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # rMLST Assembled
        logging.info('Validating rmlst_assembled.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'rmlst_assembled.csv'),
                             test_report=os.path.join(self.test_folder, 'rmlst_assembled.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)
        # serosippr
        logging.info('Validating serosippr.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'serosippr.csv'),
                             test_report=os.path.join(self.test_folder, 'serosippr.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # sistr
        logging.info('Validating sistr.tsv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'sistr.tsv'),
                             test_report=os.path.join(self.test_folder, 'sistr.tsv'),
                             columns_to_exclude=['genome'],
                             identifying_column='genome',
                             separator='\t')
        # 16S
        logging.info('Validating sixteens_full.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'sixteens_full.csv'),
                             test_report=os.path.join(self.test_folder, 'sixteens_full.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain')
        # Univec
        logging.info('Validating univec.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'univec.csv'),
                             test_report=os.path.join(self.test_folder, 'univec.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             resfinder=True)
        # Virulence
        logging.info('Validating virulence.csv')
        self.validate_report(reference_report=os.path.join(self.reference_folder, 'virulence.csv'),
                             test_report=os.path.join(self.test_folder, 'virulence.csv'),
                             columns_to_exclude=['Strain'],
                             identifying_column='Strain',
                             one_to_one=True)

    def validate_combined_metadata(self, reference_report, test_report, assembly_typer=False):
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
            self.validate_pass = False
        else:
            self.validate_pass = True

    def validate_report(self, reference_report, test_report, columns_to_exclude, identifying_column, one_to_one=False,
                        resfinder=False, separator=','):
        columns = validate.find_all_columns(csv_file=test_report,
                                            columns_to_exclude=columns_to_exclude,
                                            separator=separator)
        report_validate = validate.Validator(reference_csv=reference_report,
                                             test_csv=test_report,
                                             column_list=columns,
                                             identifying_column=identifying_column,
                                             separator=separator)
        if resfinder or one_to_one:
            if one_to_one:
                assert report_validate.check_resfinderesque_output(one_to_one=True) is True
            else:
                assert report_validate.check_resfinderesque_output() is True
        else:
            assert report_validate.all_test_columns_in_ref_and_test() is True
            assert report_validate.same_columns_in_ref_and_test() is True
            assert report_validate.check_samples_present() is True
            assert report_validate.check_columns_match() is True

    def __init__(self, reference_folder, test_folder, assembly_typer=False):
        """

        :param reference_folder:
        :param test_folder:
        :param assembly_typer:
        """
        if reference_folder.startswith('~'):
            self.reference_folder = os.path.abspath(os.path.expanduser(os.path.join(reference_folder)))
        else:
            self.reference_folder = os.path.abspath(os.path.join(reference_folder))
        assert os.path.isdir(self.reference_folder), 'Cannot locate reference report folder as specified: {rf}'\
            .format(rf=self.reference_folder)
        if test_folder.startswith('~'):
            self.test_folder = os.path.abspath(os.path.expanduser(os.path.join(test_folder)))
        else:
            self.test_folder = os.path.abspath(os.path.join(test_folder))
        assert os.path.isdir(self.test_folder), 'Cannot locate test report folder as specified: {tf}'\
            .format(tf=self.test_folder)
        self.assembly_typer = assembly_typer
        self.validate_pass = False


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Run integration tests on COWBAT pipeline')
    parser.add_argument('-r', '--reference_folder',
                        required=True,
                        help='Path to reference folder with CSV reports with expected results.')
    parser.add_argument('-t', '--test_folder',
                        required=True,
                        help='Path to test folder with CSV reports with observed results .')
    parser.add_argument('-a', '--assembly',
                        action='store_true',
                        help='The assembly typing pipeline was used to process the run, rather than full COWBAT')
    # Get the arguments into an object
    args = parser.parse_args()
    # Pretty logging!
    SetupLogging()
    # Test the reports.
    validate_outputs = ValidateCowbat(reference_folder=args.reference_folder,
                                      test_folder=args.test_folder,
                                      assembly_typer=args.assembly)
    validate_outputs.validate_cowbat()
    if validate_outputs.validate_pass:
        logging.info('COWBAT successfully validated! :D')
    else:
        logging.error('COWBAT not successfully validated.')
