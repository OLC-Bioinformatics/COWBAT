#!/usr/bin/env python
from accessoryFunctions.accessoryFunctions import SetupLogging
import os
import csv
import sys
import logging

if __name__ == '__main__':
    # Create lots of dictionaries that show what answers should be (for categorical stuff like SamplePurity or Genus)
    # or have ranges of values (for things like N50 or number of contigs - these might vary somewhat, but broadly should
    # remain the same.

    SetupLogging()
    # TODO: Add more things to check - maybe coverage depth, percent GC, AMR resistance, GeneSeekr profiles?
    # Confer with Adam/Cathy about this.

    # I see no reason for genus predictions to end up changing - these should be good to stay the same forever.
    expected_genera = {
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

    # This should also be able to stay the same, unless ConFindr gets very reworked or a new contamination detection
    # tool comes out.
    expected_contam = {
        '2013-SEQ-0132': 'Clean',
        '2014-SEQ-0136': 'Clean',
        '2014-SEQ-0276': 'Clean',
        '2014-SEQ-0933': 'Contaminated',
        '2014-SEQ-1049': 'Contaminated',
        '2015-SEQ-0423': 'Clean',
        '2015-SEQ-0626': 'Clean',
        '2015-SEQ-1361': 'Clean',
        '2017-HCLON-0380': 'Clean',
        '2017-SEQ-0905': 'Clean',
        '2017-SEQ-1501': 'Clean',
        '2018-STH-0076': 'Clean'
    }

    # I don't see any reason for these to change, barring us getting way better at somehow pulling rMLST out of
    # extremely low coverage assemblies
    expected_rmlst = {
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

    expected_sistr_cgmlst = {
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

    # These N50 ranges are very generous - unless a new assembler that does way better than SKESA/SPAdes comes out,
    # I can't see a reason that they would need to be changed.
    n50_ranges = {
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

    # These are in the same boat as N50 ranges - barring a fairly revolutionary new assembles, should be broad enough
    # that they won't need to be changed.
    contig_ranges = {
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
    total_length_ranges = {
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

    # These I can't see ever needing to change - we should have a very good idea of what depth things are.
    # Should raise a fairly major warning if this value is outside the expected ranges.
    mean_insert_size_ranges = {
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
        '2018-STH-0076': [350, 400]
    }

    # These I can't see ever needing to change - we should have a very good idea of what depth things are.
    # Should raise a fairly major warning if this value is outside the expected ranges.
    avg_cov_depth_ranges = {
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
        '2018-STH-0076': [1, 10]
    }

    # Check that correct marker genes are present.
    expected_marker_genes = {
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

    expected_vtyper_profile = {
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

    if len(sys.argv) != 2:
        logging.info('USAGE: python validate_data.py /path/to/run/folder')
        quit()

    logging.info('Validating COWBAT...')
    run_folder = sys.argv[1]

    # Start by checking that the run folder the user specified actually exists, and then check that combinedMetadata
    # exists.
    if not os.path.isdir(run_folder):
        logging.info('The run folder you specified does not exist! Verify the path to the run folder and try again.')
        quit()

    if not os.path.isfile(os.path.join(run_folder, 'reports', 'combinedMetadata.csv')):
        logging.info('Could not find combinedMetadata for run {run_folder}. '
                     'COWBAT may not have run successfully.'.format(run_folder=run_folder))
        quit()

    # Set a flag. If we make it through with no errors (flag is false) logging.info out a nice happy message.
    warning_flag = False

    # Once we know that the combinedMetadata file actually exists, time to parse through it to make sure that everything
    # looks the way it should.
    combined_metadata = os.path.join(run_folder, 'reports', 'combinedMetadata.csv')
    with open(combined_metadata) as csvfile:
        reader = csv.DictReader(csvfile)
        # For each row, verify all the things!
        for row in reader:
            # Check that items in GeneSeekr_Profile are all present.
            geneseekr_items = row['GeneSeekr_Profile'].split(';')
            for gene in expected_marker_genes[row['SeqID']]:
                if gene not in geneseekr_items:
                    logging.warning('WARNING: Expected to find {gene} in {sample}, '
                                    'but {gene} was not found in GeneSeekr_Profile'
                                    .format(gene=gene,
                                            sample=row['SeqID']))
                    warning_flag = True
            # Categorical verification: easy peasy.
            if expected_genera[row['SeqID']] != row['Genus']:
                logging.warning('WARNING: Expected genus to be {expected_genus} for {sample}, '
                                'but got genus {found_genus}!'
                                .format(expected_genus=expected_genera[row['SeqID']],
                                        sample=row['SeqID'],
                                        found_genus=row['Genus']))
                warning_flag = True
            if expected_contam[row['SeqID']] != row['SamplePurity']:
                logging.warning('WARNING: Expected SamplePurity to be {expected_purity} for {sample}, '
                                'but got purity {found_purity}!'
                                .format(expected_purity=expected_contam[row['SeqID']],
                                        sample=row['SeqID'],
                                        found_purity=row['SamplePurity']))
                warning_flag = True
            if expected_rmlst[row['SeqID']] != row['rMLST_Result']:
                logging.warning('WARNING: Expected rMLST to be {expected_rmlst} for {sample}, '
                                'but got rMLST {found_rmlst}!'
                                .format(expected_rmlst=expected_rmlst[row['SeqID']],
                                        sample=row['SeqID'],
                                        found_rmlst=row['rMLST_Result']))
                warning_flag = True
            if expected_sistr_cgmlst[row['SeqID']] != row['SISTR_serovar_cgMLST']:
                logging.warning('WARNING: Expected SISTR_serovar_cgMLST to be {expected_cgmlst} for {sample}, '
                                'but got SISTR_serovar_cgMLST {found_cgmlst}!'
                                .format(expected_cgmlst=expected_sistr_cgmlst[row['SeqID']],
                                        sample=row['SeqID'],
                                        found_cgmlst=row['SISTR_serovar_cgMLST']))
                warning_flag = True
            vtyper_items = row['Legacy_Vtyper_Profile'].split(';')
            for gene in expected_vtyper_profile[row['SeqID']]:
                if gene not in vtyper_items:
                    logging.warning('WARNING: Expected Legacy_Vtyper_Profile to be {expected_vtyper} for {sample}, '
                                    'but got Legacy_Vtyper_Profile {found_vtyper}!'
                                    .format(expected_vtyper=expected_vtyper_profile[row['SeqID']],
                                            sample=row['SeqID'],
                                            found_vtyper=row['Legacy_Vtyper_Profile']))
                    warning_flag = True
            # Now verify data that we allow ranges for (and also sometimes shows up as 0 or ND, depending on the column)
            n50 = row['N50']
            # If n50 can't be calculated because an data quality was too low for any assembly to be produced, it ends up
            # showing up as ND. If that happens, check that that's what was supposed to happen.
            if n50 == 'ND':
                if n50_ranges[row['SeqID']][0] != 'ND':
                    logging.warning('WARNING: Expected N50 did not match for {sample}. Found ND, but expected'
                                    ' {expected_n50}'
                                    .format(sample=row['SeqID'],
                                            expected_n50=n50_ranges[row['SeqID']][0]))
                    warning_flag = True
            else:
                n50 = int(n50)
                if not n50_ranges[row['SeqID']][0] <= n50 <= n50_ranges[row['SeqID']][1]:
                    logging.warning('WARNING: N50 for {sample} did not fall in expected range. N50 was {n50}. '
                                    'Expected range is from {range_bottom} to {range_top}'
                                    .format(sample=row['SeqID'],
                                            n50=n50,
                                            range_bottom=n50_ranges[row['SeqID']][0],
                                            range_top=n50_ranges[row['SeqID']][1]))
                    warning_flag = True
            # Num contigs and total length both get set to 0 instead of ND, so don't need to worry about that
            num_contigs = int(row['NumContigs'])
            if not contig_ranges[row['SeqID']][0] <= num_contigs <= contig_ranges[row['SeqID']][1]:
                logging.warning('WARNING: Number of contigs for {sample} did not fall in expected range. '
                                'Number of contigs was {num_contigs}. Expected range is from '
                                '{range_bottom} to {range_top}'
                                .format(sample=row['SeqID'],
                                        num_contigs=num_contigs,
                                        range_bottom=contig_ranges[row['SeqID']][0],
                                        range_top=contig_ranges[row['SeqID']][1]))
                warning_flag = True

            total_length = int(row['TotalLength'])
            if not total_length_ranges[row['SeqID']][0] <= total_length <= total_length_ranges[row['SeqID']][1]:
                logging.warning('WARNING: Total length for {sample} did not fall in expected range. '
                                'Total length was {total_length}. Expected range is from {range_bottom} to '
                                '{range_top}'.format(sample=row['SeqID'],
                                                     total_length=total_length,
                                                     range_bottom=total_length_ranges[row['SeqID']][0],
                                                     range_top=total_length_ranges[row['SeqID']][1]))
                warning_flag = True

            mean_ins = float(row['MeanInsertSize'])
            if not mean_insert_size_ranges[row['SeqID']][0] <= mean_ins <= mean_insert_size_ranges[row['SeqID']][1]:
                logging.warning('WARNING: Mean insert size for {sample} did not fall in expected range. '
                                'Mean insert size was {mean_insert_size}. Expected range is from {range_bottom} to '
                                '{range_top}'.format(sample=row['SeqID'],
                                                     mean_insert_size=mean_ins,
                                                     range_bottom=mean_insert_size_ranges[row['SeqID']][0],
                                                     range_top=mean_insert_size_ranges[row['SeqID']][1]))
                warning_flag = True

            avg_cov_depth = float(row['AverageCoverageDepth'])
            if not avg_cov_depth_ranges[row['SeqID']][0] <= avg_cov_depth <= avg_cov_depth_ranges[row['SeqID']][1]:
                logging.warning('WARNING: Average coverage depth for {sample} did not fall in expected range. '
                                'Average coverage depth was {avg_cov_depth}. Expected range is from {range_bottom} to '
                                '{range_top}'.format(sample=row['SeqID'],
                                                     avg_cov_depth=avg_cov_depth,
                                                     range_bottom=avg_cov_depth_ranges[row['SeqID']][0],
                                                     range_top=avg_cov_depth_ranges[row['SeqID']][1]))
                warning_flag = True
    # Also make sure that each row has the correct number of columns - sometimes in the past things have been shifted
    # over due to formatting errors of some sort.
    with open(combined_metadata) as f:
        lines = f.readlines()
    correct_column_number = len(lines[0].split(','))
    for i in range(1, len(lines)):
        if len(lines[i].split(',')) != correct_column_number:
            logging.warning('WARNING: {sample} has {bad_column_number} columns when it should have '
                            '{correct_column_number}. Check its formatting, because something has gone wrong!'
                            .format(sample=lines[i].split(',')[0],
                                    bad_column_number=len(lines[i].split(',')),
                                    correct_column_number=correct_column_number))
            warning_flag = True
    # Now that we've reached the end of our checks, let the user know if COWBAT has been successfully validated.
    if warning_flag:
        logging.warning('All checks complete, one or more warning(s) encountered. COWBAT not successfully validated :(')
    else:
        logging.info('All checks complete, no warnings encountered. COWBAT successfully validated! :D')
