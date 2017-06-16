#!/usr/bin/env python
from accessoryFunctions import *
from glob import glob
__author__ = 'adamkoziol'


class Aggregate(object):

    def handler(self):
        """
        Run all the necessary modules in the correct order
        """
        # Find all folders of the assembled run, and the names of the reports
        self.runmetadata.samples = self.reports()
        # Read all the reports, and aggregate each type into a combined report
        self.aggregate()

    def reports(self):
        """
        Find all the runs in the path. Create an object to store the metadata for each run: run name, path, name of
        reports, etc.
        :return: a list of metadata objects
        """
        printtime('Finding reports', self.start)
        # Initialise a list of runs and a list of all the metadata objects
        runs = list()
        samples = list()
        # Glob all the assembly folders into a list
        folders = glob('{}*/'.format(self.path))
        # Iterate through all the folders
        for folder in folders:
            # Ignore any previously processed aggregate reports folder
            if 'reports' not in folder:
                # Treat the external labs differently - instead of having all the assemblies in the path, each lab
                # has its own subfolder with assemblies
                if self.external:
                    subfolder = glob('{}*/'.format(folder))
                    for subsubfolder in subfolder:
                        runs.append(subsubfolder)
                else:
                    runs.append(folder)
        # Create metadata objects for each assembly run
        for assembly in sorted(runs):
            # Initialise the Metadata object
            metadata = MetadataObject()
            metadata.name = assembly.split('/')[-2] if assembly.split('/')[-3] == 'WGSspades' else '{}-{}'\
                .format(assembly.split('/')[-3], assembly.split('/')[-2])
            # Strip of any underscores and anything after an underscore
            metadata.name = metadata.name.split('_')[0]
            # Initialise Genobjects to store nested dictionaries of metadata in each metadata object
            metadata.general = GenObject()
            # Populate the Genobject
            metadata.general.path = assembly
            metadata.general.reportpath = str(assembly) + 'reports/'
            metadata.general.reports = glob('{}*.csv'.format(metadata.general.reportpath))
            # SISTR outputs a .tsv file - add it to the list of reports
            metadata.general.reports = metadata.general.reports + glob('{}*.tsv'.format(metadata.general.reportpath))
            metadata.general.reportnames = [os.path.basename(x) for x in metadata.general.reports]
            # Add all the names of the reports to the set of reports
            for reportname in metadata.general.reportnames:

                self.reportset.add(reportname)
            # Add the metadata to list of metadata objects
            samples.append(metadata)
        # Return the list of metadata objects
        return samples

    def aggregate(self):
        """
        Aggregate all reports of the same type into a master report
        """
        for report in self.reportset:
            printtime('Processing {}'.format(report.split('.')[0]), self.start)
            # Initialise the header for each report - MLST is different, as the header is different for each
            # MLST scheme. This provides a generic header instead
            header = '' if report != 'mlst.csv' else 'Strain,Genus,SequenceType,Matches,1,2,3,4,5,6,7\n'
            # Initialise a string to hold the data for each report
            data = ''
            # Open the aggregated report
            with open('{}{}'.format(self.reportpath, report), 'wb') as aggregate:
                for sample in self.runmetadata.samples:
                    # Try to open the report for this run
                    try:
                        #
                        with open('{}{}'.format(sample.general.reportpath, report), 'rb') as runreport:
                            # Only get the header from the first file
                            if not header:
                                header = runreport.readline()
                            else:
                                for row in runreport:
                                    # The final entry in a report does not have a newline character. Add \n as required
                                    if not row.endswith('\n'):
                                        row += '\n'
                                    # For certain reports, the header row is printed above each strain - ignore multiple
                                    # instances of the header
                                    if row.split(',')[0] != header.split(',')[0]:
                                        # Add the row to the string of data
                                        data += row
                    except IOError:
                        pass
                # Write the strings to the aggregate report file
                aggregate.write(header)
                aggregate.write(data)

    def __init__(self, args, starttime):
        """
        :param args: Command line arguments
        :param starttime: time the script was started
        """
        self.args = args
        self.start = starttime
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.external = args.external
        # Initialise class variables
        self.samples = list()
        self.reportset = set()
        self.runmetadata = MetadataObject()
        self.reportpath = os.path.join(self.path, 'reports', '')
        # Make the report folder
        make_path(self.reportpath)
        # Run the analyses
        self.handler()

if __name__ == '__main__':
    import time
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Combine all the reports created by the pipeline over several runs into a '
                                        'single report for each type')
    parser.add_argument('path',
                        help='Input directory. REQUIRED')
    parser.add_argument('-e', '--external',
                        action='store_true',
                        help='The directory format for local analyses is different than for the external analyses. '
                             'Add this flag if the external analyses are to be processed')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Define the start time
    start = time.time()
    # Run the script
    Aggregate(arguments, start)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')
