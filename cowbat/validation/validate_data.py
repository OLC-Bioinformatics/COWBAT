#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import SetupLogging
from argparse import ArgumentParser
import logging
import csv
import os


class Validate(object):

    def main(self):
        self.validation()
        self.validate_columns()
        self.warnings()

    def validation(self):
        """

        :return:
        """
        logging.info('Validating COWBAT.')
        # Open the combined metadata report
        with open(self.combined_metadata) as csvfile:
            reader = csv.DictReader(csvfile)
            # For each row, verify all the things!
            for row in reader:
                # Check that items in GeneSeekr_Profile are all present.
                geneseekr_items = row['GeneSeekr_Profile'].split(';')
                for gene in self.expected_marker_genes[row['SeqID']]:
                    if gene not in geneseekr_items:
                        logging.warning('WARNING: Expected to find {gene} in {sample}, '
                                        'but {gene} was not found in GeneSeekr_Profile'
                                        .format(gene=gene,
                                                sample=row['SeqID']))
                        self.warning_flag = True
                # Categorical verification
                # Expected genera
                if self.expected_genera[row['SeqID']] != row['Genus']:
                    logging.warning('WARNING: Expected genus to be {expected_genus} for {sample}, '
                                    'but got genus {found_genus}!'
                                    .format(expected_genus=self.expected_genera[row['SeqID']],
                                            sample=row['SeqID'],
                                            found_genus=row['Genus']))
                    self.warning_flag = True
                # ConFindr
                if self.expected_contam[row['SeqID']] != row['SamplePurity']:
                    logging.warning('WARNING: Expected SamplePurity to be {expected_purity} for {sample}, '
                                    'but got purity {found_purity}!'
                                    .format(expected_purity=self.expected_contam[row['SeqID']],
                                            sample=row['SeqID'],
                                            found_purity=row['SamplePurity']))
                    self.warning_flag = True
                # rMLST
                if self.expected_rmlst[row['SeqID']] != row['rMLST_Result']:
                    logging.warning('WARNING: Expected rMLST to be {expected_rmlst} for {sample}, '
                                    'but got rMLST {found_rmlst}!'
                                    .format(expected_rmlst=self.expected_rmlst[row['SeqID']],
                                            sample=row['SeqID'],
                                            found_rmlst=row['rMLST_Result']))
                    self.warning_flag = True
                # SISRE cgMLST
                if self.expected_sistr_cgmlst[row['SeqID']] != row['SISTR_serovar_cgMLST']:
                    logging.warning('WARNING: Expected SISTR_serovar_cgMLST to be {expected_cgmlst} for {sample}, '
                                    'but got SISTR_serovar_cgMLST {found_cgmlst}!'
                                    .format(expected_cgmlst=self.expected_sistr_cgmlst[row['SeqID']],
                                            sample=row['SeqID'],
                                            found_cgmlst=row['SISTR_serovar_cgMLST']))
                    self.warning_flag = True
                # Vtyper
                vtyper_items = row['Vtyper_Profile'].split(';')
                for gene in self.expected_vtyper_profile[row['SeqID']]:
                    if gene not in vtyper_items:
                        logging.warning('WARNING: Expected Legacy_Vtyper_Profile to be {expected_vtyper} for {sample}, '
                                        'but got Legacy_Vtyper_Profile {found_vtyper}!'
                                        .format(expected_vtyper=self.expected_vtyper_profile[row['SeqID']],
                                                sample=row['SeqID'],
                                                found_vtyper=row['Legacy_Vtyper_Profile']))
                        self.warning_flag = True

                # Verify data that we allow ranges for (and also sometimes shows up as 0 or ND, depending on the column)
                n50 = row['N50']
                # If n50 can't be calculated because an data quality was too low for any assembly to be produced,
                # it ends up showing up as ND. If that happens, check that that's what was supposed to happen.
                if n50 == 'ND':
                    if self.n50_ranges[row['SeqID']][0] != 'ND' and self.n50_ranges[row['SeqID']][0] != 0:
                        logging.warning('WARNING: Expected N50 did not match for {sample}. Found ND, but expected'
                                        ' {expected_n50}'
                                        .format(sample=row['SeqID'],
                                                expected_n50=self.n50_ranges[row['SeqID']][0]))
                        self.warning_flag = True
                else:
                    n50 = int(n50)
                    if not self.n50_ranges[row['SeqID']][0] <= n50 <= self.n50_ranges[row['SeqID']][1]:
                        logging.warning('WARNING: N50 for {sample} did not fall in expected range. N50 was {n50}. '
                                        'Expected range is from {range_bottom} to {range_top}'
                                        .format(sample=row['SeqID'],
                                                n50=n50,
                                                range_bottom=self.n50_ranges[row['SeqID']][0],
                                                range_top=self.n50_ranges[row['SeqID']][1]))
                        self.warning_flag = True
                # Num contigs and total length both get set to 0 instead of ND, so don't need to worry about that
                num_contigs = int(row['NumContigs'])
                if not self.contig_ranges[row['SeqID']][0] <= num_contigs <= self.contig_ranges[row['SeqID']][1]:
                    logging.warning('WARNING: Number of contigs for {sample} did not fall in expected range. '
                                    'Number of contigs was {num_contigs}. Expected range is from '
                                    '{range_bottom} to {range_top}'
                                    .format(sample=row['SeqID'],
                                            num_contigs=num_contigs,
                                            range_bottom=self.contig_ranges[row['SeqID']][0],
                                            range_top=self.contig_ranges[row['SeqID']][1]))
                    self.warning_flag = True
                # Total length
                total_length = int(row['TotalLength'])
                if not self.total_length_ranges[row['SeqID']][0] <= total_length <= self.total_length_ranges[row['SeqID']][1]:
                    logging.warning('WARNING: Total length for {sample} did not fall in expected range. '
                                    'Total length was {total_length}. Expected range is from {range_bottom} to '
                                    '{range_top}'.format(sample=row['SeqID'],
                                                         total_length=total_length,
                                                         range_bottom=self.total_length_ranges[row['SeqID']][0],
                                                         range_top=self.total_length_ranges[row['SeqID']][1]))
                    self.warning_flag = True
                # Mean insert size
                try:
                    mean_ins = float(row['MeanInsertSize'])
                except ValueError:
                    mean_ins = 0
                if not self.mean_insert_size_ranges[row['SeqID']][0] <= mean_ins <= self.mean_insert_size_ranges[row['SeqID']][1]:
                    logging.warning('WARNING: Mean insert size for {sample} did not fall in expected range. '
                                    'Mean insert size was {mean_insert_size}. Expected range is from {range_bottom} to '
                                    '{range_top}'.format(sample=row['SeqID'],
                                                         mean_insert_size=mean_ins,
                                                         range_bottom=self.mean_insert_size_ranges[row['SeqID']][0],
                                                         range_top=self.mean_insert_size_ranges[row['SeqID']][1]))
                    self.warning_flag = True
                # Average coverage depth
                try:
                    avg_cov_depth = float(row['AverageCoverageDepth'])
                except ValueError:
                    avg_cov_depth = 0
                if not self.avg_cov_depth_ranges[row['SeqID']][0] <= avg_cov_depth <= self.avg_cov_depth_ranges[row['SeqID']][1]:
                    logging.warning('WARNING: Average coverage depth for {sample} did not fall in expected range. '
                                    'Average coverage depth was {avg_cov_depth}. Expected range is from {range_b} to '
                                    '{range_t}'.format(sample=row['SeqID'],
                                                       avg_cov_depth=avg_cov_depth,
                                                       range_b=self.avg_cov_depth_ranges[row['SeqID']][0],
                                                       range_t=self.avg_cov_depth_ranges[row['SeqID']][1]))
                    self.warning_flag = True
                # Number of predicted genes
                num_predicted_genes = float(row['TotalPredictedGenes'])
                if not self.num_predicted_genes[row['SeqID']][0] <= num_predicted_genes <= \
                       self.num_predicted_genes[row['SeqID']][1]:
                    logging.warning('WARNING: Number of predicted genes for {sample} did not fall in expected range. '
                                    'Number of predicted genes was {num_genes}. Expected range is from {range_b} to '
                                    '{range_t}'.format(sample=row['SeqID'],
                                                       num_genes=num_predicted_genes,
                                                       range_b=self.num_predicted_genes[row['SeqID']][0],
                                                       range_t=self.num_predicted_genes[row['SeqID']][1]))
                    self.warning_flag = True
                # Number MASH matching hashes
                try:
                    matching_hashes = float(row['MASH_NumMatchingHashes'].split('/')[0])
                except ValueError:
                    matching_hashes = 0
                if not self.mash_matching_hashes[row['SeqID']][0] <= matching_hashes <= \
                       self.mash_matching_hashes[row['SeqID']][1]:
                    logging.warning(
                        'WARNING: Number of MASH matching hashes for {sample} did not fall in expected range. '
                        'Number of MASH matching hashes was {num_hashes}. Expected range is from {range_b} to '
                        '{range_t}'.format(sample=row['SeqID'],
                                           num_hashes=matching_hashes,
                                           range_b=self.mash_matching_hashes[row['SeqID']][0],
                                           range_t=self.mash_matching_hashes[row['SeqID']][1]))
                    self.warning_flag = True

    def validate_columns(self):
        """
        Ensure that each row has the correct number of columns. Sometimes in the past things have been shifted
        over due to formatting errors of some sort.
        """
        with open(self.combined_metadata) as f:
            lines = f.readlines()
        correct_column_number = len(lines[0].split(','))
        for i in range(1, len(lines)):
            if len(lines[i].split(',')) != correct_column_number:
                logging.warning('WARNING: {sample} has {bad_column_number} columns when it should have '
                                '{correct_column_number}. Check its formatting, because something has gone wrong!'
                                .format(sample=lines[i].split(',')[0],
                                        bad_column_number=len(lines[i].split(',')),
                                        correct_column_number=correct_column_number))
                self.warning_flag = True

    def warnings(self):
        """
        Now that we've reached the end of our checks, let the user know if COWBAT has been successfully validated.
        """
        if self.warning_flag:
            logging.warning('All checks complete, one or more warning(s) encountered. '
                            'COWBAT not successfully validated :(')
        else:
            logging.info('All checks complete, no warnings encountered. COWBAT successfully validated! :D')

    def __init__(self, run_folder, assembly_typing):
        """

        :param run_folder:
        """
        SetupLogging()
        self.run_folder = os.path.join(run_folder)
        assert os.path.isdir(self.run_folder), 'Cannot locate the specified report path: {rp}'\
            .format(rp=self.run_folder)
        self.combined_metadata = os.path.join(self.run_folder, 'reports', 'combinedMetadata.csv')
        # Set a flag. If we make it through with no errors (flag is false) logging.info out a nice happy message.
        self.warning_flag = False
        self.assembly = assembly_typing
        # Create dictionaries that show what answers should be (for categorical stuff like SamplePurity or Genus) or
        # have ranges of values (for things like N50 or number of contigs. These might vary somewhat, but broadly should
        # remain the same.
        # TODO: Add more things to check - maybe coverage depth, percent GC, AMR resistance, GeneSeekr profiles?
        # Confer with Adam/Cathy about this.

        # I see no reason for genus predictions to end up changing - these should be good to stay the same forever.
        self.expected_genera = {
            '2013-SEQ-0132': 'Escherichia',
            '2014-SEQ-0136': 'Listeria',
            '2014-SEQ-0276': 'Escherichia',
            '2014-SEQ-0933': 'Listeria',
            '2014-SEQ-1049': 'Salmonella',
            '2015-SEQ-0423': 'Salmonella',
            '2015-SEQ-0626': 'Pseudomonas',
            '2015-SEQ-1361': 'Moellerella',
            '2017-HCLON-0380': 'Escherichia',
            '2017-SEQ-0905': 'Acinetobacter',
            '2017-SEQ-1501': 'Listeria',
            '2018-STH-0076': 'ND'
        }
        # I don't see any reason for these to change, barring us getting way better at somehow pulling rMLST out of
        # extremely low coverage assemblies
        self.expected_rmlst = {
            '2013-SEQ-0132': 'new',
            '2014-SEQ-0136': '16022',
            '2014-SEQ-0276': '2124',
            '2014-SEQ-0933': '15987',
            '2014-SEQ-1049': '3811',
            '2015-SEQ-0423': 'new',
            '2015-SEQ-0626': 'new',
            '2015-SEQ-1361': '39701',
            '2017-HCLON-0380': '2011',
            '2017-SEQ-0905': '9659',
            '2017-SEQ-1501': '32448',
            '2018-STH-0076': 'new'
        }
        self.expected_sistr_cgmlst = {
            '2013-SEQ-0132': 'ND',
            '2014-SEQ-0136': 'ND',
            '2014-SEQ-0276': 'ND',
            '2014-SEQ-0933': 'ND',
            '2014-SEQ-1049': 'Berta',
            '2015-SEQ-0423': 'Infantis',
            '2015-SEQ-0626': 'ND',
            '2015-SEQ-1361': 'ND',
            '2017-HCLON-0380': 'ND',
            '2017-SEQ-0905': 'ND',
            '2017-SEQ-1501': 'ND',
            '2018-STH-0076': 'ND'
        }
        self.expected_marker_genes = {
            '2013-SEQ-0132': ['uidA', 'VT2'],
            '2014-SEQ-0136': ['hlyALm', 'IGS', 'inlJ'],
            '2014-SEQ-0276': ['eae', 'O157', 'VT1', 'VT2', 'uidA'],
            '2014-SEQ-0933': ['hlyALm', 'IGS', 'inlJ'],
            '2014-SEQ-1049': ['invA', 'stn'],
            '2015-SEQ-0423': ['invA', 'stn'],
            '2015-SEQ-0626': ['ND'],
            '2015-SEQ-1361': ['ND'],
            '2017-HCLON-0380': ['uidA'],
            '2017-SEQ-0905': ['ND'],
            '2017-SEQ-1501': ['hlyALm', 'IGS', 'inlJ'],
            '2018-STH-0076': ['ND']
        }
        self.expected_vtyper_profile = {
            '2013-SEQ-0132': ['vtx2a'],
            '2014-SEQ-0136': ['ND'],
            '2014-SEQ-0276': ['vtx1a', 'vtx2a'],
            '2014-SEQ-0933': ['ND'],
            '2014-SEQ-1049': ['ND'],
            '2015-SEQ-0423': ['ND'],
            '2015-SEQ-0626': ['ND'],
            '2015-SEQ-1361': ['ND'],
            '2017-HCLON-0380': ['ND'],
            '2017-SEQ-0905': ['ND'],
            '2017-SEQ-1501': ['ND'],
            '2018-STH-0076': ['ND']
        }
        # These N50 ranges are very generous - unless a new assembler that does way better than SKESA/SPAdes comes out,
        # I can't see a reason that they would need to be changed.
        self.n50_ranges = {
            '2013-SEQ-0132': [1000, 5000],
            '2014-SEQ-0136': [15000, 30000],
            '2014-SEQ-0276': [100000, 200000],
            '2014-SEQ-0933': [400000, 600000],
            '2014-SEQ-1049': [125000, 225000],
            '2015-SEQ-0423': [1000, 5000],
            '2015-SEQ-0626': [150000, 250000],
            '2015-SEQ-1361': [150000, 250000],
            '2017-HCLON-0380': [80000, 140000],
            '2017-SEQ-0905': [30000, 60000],
            '2017-SEQ-1501': [200000, 375000],
            '2018-STH-0076': [0, 500]
        }
        # These are in the same boat as N50 ranges. Barring a fairly revolutionary new assembler, should be broad enough
        # that they won't need to be changed.
        self.contig_ranges = {
            '2013-SEQ-0132': [1000, 5000],
            '2014-SEQ-0136': [190, 250],
            '2014-SEQ-0276': [100, 200],
            '2014-SEQ-0933': [10, 30],
            '2014-SEQ-1049': [40, 80],
            '2015-SEQ-0423': [1000, 5000],
            '2015-SEQ-0626': [80, 140],
            '2015-SEQ-1361': [30, 70],
            '2017-HCLON-0380': [150, 325],
            '2017-SEQ-0905': [140, 230],
            '2017-SEQ-1501': [15, 45],
            '2018-STH-0076': [0, 50]
        }
        # These I can't see ever needing to change - we should have a very good idea of what length things are.
        # Should raise a fairly major warning if range limits are ever exceeded.
        self.total_length_ranges = {
            '2013-SEQ-0132': [4600000, 5000000],
            '2014-SEQ-0136': [2800000, 3000000],
            '2014-SEQ-0276': [5100000, 5500000],
            '2014-SEQ-0933': [2900000, 3100000],
            '2014-SEQ-1049': [4600000, 4900000],
            '2015-SEQ-0423': [3900000, 4200000],
            '2015-SEQ-0626': [6100000, 6600000],
            '2015-SEQ-1361': [3100000, 3600000],
            '2017-HCLON-0380': [5100000, 5500000],
            '2017-SEQ-0905': [3700000, 4100000],
            '2017-SEQ-1501': [2900000, 3100000],
            '2018-STH-0076': [0, 10000]
        }
        self.num_predicted_genes = {
            '2013-SEQ-0132': [6200, 6400],
            '2014-SEQ-0136': [3000, 3200],
            '2014-SEQ-0276': [5100, 5300],
            '2014-SEQ-0933': [2950, 3150],
            '2014-SEQ-1049': [4350, 4550],
            '2015-SEQ-0423': [5700, 5900],
            '2015-SEQ-0626': [5700, 5900],
            '2015-SEQ-1361': [2900, 3100],
            '2017-HCLON-0380': [5050, 5250],
            '2017-SEQ-0905': [3600, 3800],
            '2017-SEQ-1501': [2850, 3050],
            '2018-STH-0076': [0, 100]
        }
        # These dictionaries look different depending on whether the assembly typing pipeline or COWBAT was used
        if not self.assembly:
            # This should also be able to stay the same, unless ConFindr gets very reworked or a new contamination
            # detection tool comes out.
            self.expected_contam = {
                '2013-SEQ-0132': 'Clean',
                '2014-SEQ-0136': 'Clean',
                '2014-SEQ-0276': 'Clean',
                '2014-SEQ-0933': 'Clean',
                '2014-SEQ-1049': 'Clean',
                '2015-SEQ-0423': 'Clean',
                '2015-SEQ-0626': 'Contaminated',
                '2015-SEQ-1361': 'Clean',
                '2017-HCLON-0380': 'Clean',
                '2017-SEQ-0905': 'Clean',
                '2017-SEQ-1501': 'Clean',
                '2018-STH-0076': 'Clean'
            }
            # These I can't see ever needing to change - we should have a very good idea of what depth things are.
            # Should raise a fairly major warning if this value is outside the expected ranges.
            self.mean_insert_size_ranges = {
                '2013-SEQ-0132': [400, 450],
                '2014-SEQ-0136': [325, 375],
                '2014-SEQ-0276': [550, 600],
                '2014-SEQ-0933': [300, 350],
                '2014-SEQ-1049': [425, 475],
                '2015-SEQ-0423': [375, 425],
                '2015-SEQ-0626': [500, 550],
                '2015-SEQ-1361': [400, 450],
                '2017-HCLON-0380': [350, 400],
                '2017-SEQ-0905': [250, 275],
                '2017-SEQ-1501': [375, 425],
                '2018-STH-0076': [0, 400]
            }
            # These I can't see ever needing to change - we should have a very good idea of what depth things are.
            # Should raise a fairly major warning if this value is outside the expected ranges.
            self.avg_cov_depth_ranges = {
                '2013-SEQ-0132': [5, 20],
                '2014-SEQ-0136': [30, 50],
                '2014-SEQ-0276': [90, 110],
                '2014-SEQ-0933': [140, 160],
                '2014-SEQ-1049': [120, 140],
                '2015-SEQ-0423': [1, 20],
                '2015-SEQ-0626': [70, 90],
                '2015-SEQ-1361': [140, 160],
                '2017-HCLON-0380': [55, 75],
                '2017-SEQ-0905': [95, 115],
                '2017-SEQ-1501': [125, 145],
                '2018-STH-0076': [0, 10]
            }
            self.mash_matching_hashes = {
                '2013-SEQ-0132': [900, 950],
                '2014-SEQ-0136': [900, 950],
                '2014-SEQ-0276': [835, 885],
                '2014-SEQ-0933': [865, 915],
                '2014-SEQ-1049': [900, 950],
                '2015-SEQ-0423': [865, 915],
                '2015-SEQ-0626': [350, 400],
                '2015-SEQ-1361': [800, 850],
                '2017-HCLON-0380': [675, 725],
                '2017-SEQ-0905': [800, 850],
                '2017-SEQ-1501': [900, 950],
                '2018-STH-0076': [0, 50]
            }
        else:
            self.expected_contam = {
                '2013-SEQ-0132': 'ND',
                '2014-SEQ-0136': 'ND',
                '2014-SEQ-0276': 'ND',
                '2014-SEQ-0933': 'ND',
                '2014-SEQ-1049': 'ND',
                '2015-SEQ-0423': 'ND',
                '2015-SEQ-0626': 'ND',
                '2015-SEQ-1361': 'ND',
                '2017-HCLON-0380': 'ND',
                '2017-SEQ-0905': 'ND',
                '2017-SEQ-1501': 'ND',
                '2018-STH-0076': 'ND'
            }
            # These I can't see ever needing to change - we should have a very good idea of what depth things are.
            # Should raise a fairly major warning if this value is outside the expected ranges.
            self.mean_insert_size_ranges = {
                '2013-SEQ-0132': [0, 0],
                '2014-SEQ-0136': [0, 0],
                '2014-SEQ-0276': [0, 0],
                '2014-SEQ-0933': [0, 0],
                '2014-SEQ-1049': [0, 0],
                '2015-SEQ-0423': [0, 0],
                '2015-SEQ-0626': [0, 0],
                '2015-SEQ-1361': [0, 0],
                '2017-HCLON-0380': [0, 0],
                '2017-SEQ-0905': [0, 0],
                '2017-SEQ-1501': [0, 0],
                '2018-STH-0076': [0, 0]
            }
            # These I can't see ever needing to change - we should have a very good idea of what depth things are.
            # Should raise a fairly major warning if this value is outside the expected ranges.
            self.avg_cov_depth_ranges = {
                '2013-SEQ-0132': [0, 0],
                '2014-SEQ-0136': [0, 0],
                '2014-SEQ-0276': [0, 0],
                '2014-SEQ-0933': [0, 0],
                '2014-SEQ-1049': [0, 0],
                '2015-SEQ-0423': [0, 0],
                '2015-SEQ-0626': [0, 0],
                '2015-SEQ-1361': [0, 0],
                '2017-HCLON-0380': [0, 0],
                '2017-SEQ-0905': [0, 0],
                '2017-SEQ-1501': [0, 0],
                '2018-STH-0076': [0, 0]
            }
            self.mash_matching_hashes = {
                '2013-SEQ-0132': [850, 950],
                '2014-SEQ-0136': [900, 1000],
                '2014-SEQ-0276': [900, 1000],
                '2014-SEQ-0933': [900, 1000],
                '2014-SEQ-1049': [900, 1000],
                '2015-SEQ-0423': [800, 905],
                '2015-SEQ-0626': [450, 550],
                '2015-SEQ-1361': [900, 1000],
                '2017-HCLON-0380': [725, 825],
                '2017-SEQ-0905': [825, 925],
                '2017-SEQ-1501': [900, 1000],
                '2018-STH-0076': [0, 50]
            }


if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Run integration tests on COWBAT pipeline')
    parser.add_argument('-r', '--run_folder',
                        required=True,
                        help='Provide the location of the folder containing the pipeline reports')
    parser.add_argument('-a', '--assembly',
                        action='store_true',
                        default=False,
                        help='The assembly typing pipeline was used to process the run, rather than full COWBAT')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Run the pipeline
    validation = Validate(run_folder=arguments.run_folder,
                          assembly_typing=arguments.assembly)
    validation.main()
