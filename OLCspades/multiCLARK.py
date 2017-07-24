#!/usr/bin/env python
from accessoryFunctions import *
from glob import glob
import errno
import automateCLARK
__author__ = 'adamkoziol'


class Multi(object):

    def objectifier(self):
        """
        Create objects to be passed to automateCLARK
        """
        printtime('Finding sequencing runs', self.start)
        for lab in self.labs:
            # Find all the sequencing runs stored in each lab folder
            labfolder = os.path.join(self.path, lab, '')
            seqrun = sorted(glob('{}*/'.format(labfolder)))
            for subfolder in seqrun:
                # Initialise the Metadata object
                metadata = MetadataObject()
                metadata.name = '{}-{}'.format(subfolder.split('/')[-3], subfolder.split('/')[-2])
                # Initialise Genobjects to store nested dictionaries of metadata in each metadata object
                metadata.general = GenObject()
                # Set the path and sequence path to be passed to automateCLARK
                metadata.path = os.path.join(self.outputpath, metadata.name, '')
                metadata.sequencepath = os.path.join(metadata.path, 'sequences', '')
                # Create these paths
                make_path(metadata.sequencepath)
                # Link all the .fastq.gz files to the sequence path
                sequences = sorted(glob('{}*.gz'.format(subfolder)))
                # Creating relative symlinks requires that the links be created from within the source directory
                os.chdir(subfolder)
                for sequence in sequences:
                    try:
                        # Create the relative symlink
                        os.symlink(sequence, '{}/{}'.format(os.path.relpath(metadata.sequencepath, subfolder),
                                                            os.path.basename(sequence)))
                    # Except os errors
                    except OSError as exception:
                        # If the os error is anything but directory exists, then raise
                        if exception.errno != errno.EEXIST:
                            raise
                # Set the required attributes for automateCLARK
                metadata.databasepath = self.databasepath
                metadata.clarkpath = self.clarkpath
                metadata.rank = self.rank
                metadata.database = self.database
                metadata.cutoff = self.cutoff
                metadata.threads = self.cpus
                metadata.filter = self.filter
                # Add the metadata object to the list of objects
                self.runmetadata.append(metadata)
                printtime('Processing {}'.format(metadata.name), self.start)
                # Call automateCLARK, and pass the metadata object as well as other variables
                self.clarkdata.append(automateCLARK.CLARK(metadata, self.commit, self.start, self.homepath))

    def aggregator(self):
        """
        Combine all the reports from the individual sequencing runs
        """
        import xlsxwriter
        printtime('Creating combined report', self.start)
        from csv import DictReader
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook('{}/abundance.xlsx'.format(self.reportpath))
        make_path(self.reportpath)
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 8
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
        bold.set_align('center')
        # Format for data cells. Monotype, size 8, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
        courier.set_align('top')
        # Set the custom width for 5 and 6 to be 15
        worksheet.set_column(5, 5, 15)
        worksheet.set_column(6, 6, 20)
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        # Determine the analysis type - .fasta or .fastq
        for run in self.clarkdata:
            # All the analyses should be the same, so just use the value from the last iteration
            self.extension = run.runmetadata.extension
        # List of the headers to use
        headers = ['Strain', 'Name', 'TaxID', 'Lineage', 'Count', 'Proportion_All(%)', 'Proportion_Classified(%)']
        # Add an additional header for .fasta analyses
        if self.extension == '.fasta':
            headers.insert(4, 'TotalBP')
        # Populate the headers
        for category in headers:
            # Write the data in the specified cell (row, col) using the bold format
            worksheet.write(row, col, category, bold)
            # Move to the next column to write the next category
            col += 1
        # Data starts in row 1
        row = 1
        # Initialise variables to hold the longest names; used in setting the column width
        longeststrain = 0
        longestname = 0
        longestlineage = 0
        # Extract all the taxonomic groups that pass the cutoff from the abundance file
        for run in self.clarkdata:
            for sample in run.runmetadata.samples:
                # Every record starts at column 0
                col = 0
                # Write the strain name
                worksheet.write(row, col, sample.name, courier)
                col += 1
                # Initialise a dictionary to store the species above the cutoff in the sample
                sample.general.passfilter = list()
                # Abundance file as a dictionary
                abundancedict = DictReader(open(sample.general.abundance))
                # Filter abundance to taxIDs with at least self.cutoff% of the total proportion
                for result in abundancedict:
                    # The UNKNOWN category doesn't contain a 'Lineage' column, and therefore, subsequent columns are
                    # shifted out of proper alignment, and do not contain the appropriate data
                    try:
                        if float(result['Proportion_All(%)']) > self.cutoff:
                            sample.general.passfilter.append(result)
                    except ValueError:
                        pass
                # Determine the longest name of all the strains, and use it to set the width of column 0
                if len(sample.name) > longeststrain:
                    longeststrain = len(sample.name)
                    worksheet.set_column(0, 0, longeststrain)
                # Sort the abundance results based on the highest count
                sortedabundance = sorted(sample.general.passfilter, key=lambda x: int(x['Count']), reverse=True)
                # Set of contigs from the classification file. For some reason, certain contigs are represented multiple
                # times in the classification file. As far as I can tell, these multiple representations are always
                # classified the same, and, therefore, should be treated as duplicates, and ignored
                contigset = set()
                for result in sortedabundance:
                    # Add the total number of base pairs classified for each TaxID. As only the total number of contigs
                    # classified as a particular TaxID are presented in the report, it can be confusing many
                    # small contigs are classified to a particular TaxID e.g. 56 contigs map to TaxID 28901, and 50
                    # contigs map to TaxID 630, however, added together, those 56 contigs are 4705838 bp, while the 50
                    # contigs added together are only 69602 bp. While this is unlikely a pure culture, only
                    # 69602 / (4705838 + 69602) = 1.5% of the total bp map to TaxID 630 compared to 45% of the contigs
                    if run.runmetadata.extension == '.fasta':
                        # Initialise a variable to store the total bp mapped to the TaxID
                        result['TotalBP'] = int()
                        # Read the classification file into a dictionary
                        classificationdict = DictReader(open(sample.general.classification))
                        # Read through each contig classification in the dictionary
                        for contig in classificationdict:
                            # Pull out each contig with a TaxID that matches the TaxID of the result of interest, and
                            # is not present in a set of contigs that have already been added to the dictionary
                            if result['TaxID'] == contig[' Assignment'] and contig['Object_ID'] not in contigset:
                                # Increment the total bp mapping to the TaxID by the integer of each contig
                                result['TotalBP'] += int(contig[' Length'])
                                # Avoid duplicates by adding the contig name to the set of contigs
                                contigset.add(contig['Object_ID'])
                    # Print the results to file
                    # Ignore the first header, as it is the strain name, which has already been added to the report
                    dictionaryheaders = headers[1:]
                    for header in dictionaryheaders:
                        data = result[header]
                        worksheet.write(row, col, data, courier)
                        col += 1
                        # Determine the longest name of all the matches, and use it to set the width of column 0
                        if len(result['Name']) > longestname:
                            longestname = len(result['Name'])
                            worksheet.set_column(1, 1, longestname)
                        # Do the same for the lineages
                        if len(result['Lineage']) > longestlineage:
                            longestlineage = len(result['Lineage'])
                            worksheet.set_column(3, 3, longestlineage)
                    # Increase the row
                    row += 1
                    # Set the column to 1
                    col = 1
        # Close the workbook
        workbook.close()

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: Command line arguments
        :param startingtime: time the script was started
        """
        import multiprocessing
        self.args = args
        self.commit = str(pipelinecommit)
        self.start = startingtime
        self.homepath = scriptpath
        # Set variables from the arguments
        self.path = os.path.join(args.path, '')
        assert os.path.isdir(self.path), u'Supplied path is not a valid directory {0!r:s}'.format(self.path)
        self.outputpath = os.path.join(args.outputpath, '')
        assert os.path.isdir(self.outputpath), u'Supplied output path is not a valid directory {0!r:s}'\
            .format(self.outputpath)
        self.clarkpath = os.path.join(args.clarkpath, '')
        assert os.path.isdir(self.clarkpath), u'Supplied path to CLARK executables is not a valid directory {0!r:s}' \
            .format(self.clarkpath)
        self.databasepath = os.path.join(args.databasepath, '')
        assert os.path.isdir(self.databasepath), u'Supplied path to CLARK database is not a valid directory {0!r:s}' \
            .format(self.databasepath)
        self.database = args.database
        self.rank = args.rank
        self.cutoff = float(args.cutoff) * 100
        self.filter = args.filter
        # Create metadata objects to store data for each sequencing run
        self.runmetadata = list()
        self.clarkdata = list()
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = int(args.threads if args.threads else multiprocessing.cpu_count())
        # Create class variables
        self.reportpath = os.path.join(self.outputpath, 'reports')
        self.extension = str()
        # A list of external labs to be included in the analyses - add or subtract as required
        self.labs = ['BUR', 'CAL', 'DAR', 'GTA', 'STH']
        # Run the analyses
        self.objectifier()
        # Create a combined report
        self.aggregator()

if __name__ == '__main__':
    import time
    import subprocess
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Run CLARK analyses on all sequencing runs within a set folder hierarchy')
    parser.add_argument('path',
                        help='Input directory. REQUIRED')
    parser.add_argument('-o', '--outputpath',
                        required=True,
                        help='Path in which to store the outputs')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='Path of CLARK database files to use.')
    parser.add_argument('-C', '--clarkpath',
                        required=True,
                        help='Path to the CLARK scripts')
    parser.add_argument('-r', '--rank',
                        default='species',
                        help='Choose the taxonomic rank to use in the analysis: species (the default value), genus, '
                             'family, order, class or phylum')
    parser.add_argument('-D', '--database',
                        default='bacteria',
                        help='Choose the database to use in the analysis (one or more from: bacteria, viruses, human,'
                             'and custom. To select more than one, use commas to separate your selection. Custom'
                             'databases need to be set up before they will work')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-f', '--filter',
                        action='store_true',
                        help='Optionally split the samples based on taxonomic assignment')
    parser.add_argument('-c', '--cutoff',
                        default=0.01,
                        help='Cutoff value for setting taxIDs to use when filtering fastq files. Defaults to 1 percent.'
                             ' Please note that you must use a decimal format: enter 0.05 to get a 5 percent cutoff.')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Define the start time
    start = time.time()
    # Run the script
    Multi(arguments, commit, start, homepath)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - start) + '\033[0m')
