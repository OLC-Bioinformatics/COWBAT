#!/usr/bin/env python
from accessoryFunctions import printtime
from glob import glob
import os
__author__ = 'adamkoziol'


class Offhours(object):

    def assertpathsandfiles(self):
        """Assertions to make sure that arguments are at least mostly valid"""
        # Assertion to ensure that the MiSeq path exists
        assert os.path.isdir(self.miseqpath), u'MiSeqPath is not a valid directory {0!r:s}'.format(self.miseqpath)
        # If the miseq folder name is not provided, the default of the most recent run will be used
        if not self.miseqfolder:
            # Get a list of folders
            miseqfolders = glob('{}*/'.format(self.miseqpath))
            self.miseqfolder = sorted(miseqfolders)[-1]
            # Create :miseqfoldername to store the name of this folder by splitting the path and taking the second
            # last piece (it's not the last piece because the folder has a trailing slash)
            self.miseqfoldername = self.miseqfolder.split("/")[-2]
        # Otherwise add the folder to the miseq path to yield the destination folder
        else:
            # Set the folder name before adding the path to the miseq path
            self.miseqfoldername = self.miseqfolder
            self.miseqfolder = self.miseqpath + self.miseqfolder + "/"
            # Assert to ensure that the folder exists
            assert os.path.isdir(self.miseqfolder), u'MiSeqFolder is not a valid directory {0!r:s}'\
                .format(self.miseqfolder)
        # Pull the data from the SampleSheet.csv
        if self.customsamplesheet:
            self.samplesheet = self.customsamplesheet
            assert os.path.isfile(self.customsamplesheet), u'Could not find CustomSampleSheet as entered: {0!r:s}'\
                .format(self.customsamplesheet)
        # Otherwise use the SampleSheet.csv located in :self.miseqfolder
        else:
            self.samplesheet = self.miseqfolder + "/SampleSheet.csv"

    def numberofsamples(self):
        """Count the number of samples is the samplesheet"""
        # Initialise variables to store line data
        idline = 0
        linenumber = 0
        # Parse the sample sheet to find the number of samples
        with open(self.samplesheet, "rb") as ssheet:
            # Use enumerate to iterate through the lines in the sample sheet to retrieve the line number and the data
            for linenumber, entry in enumerate(ssheet):
                # Once Sample_ID is encountered
                if "Sample_ID" in entry:
                    # Set the id line as the current line number
                    idline = linenumber
        # :samplecount is the last line number in the file minus the line number of Sample_ID
        self.samplecount = linenumber - idline
        printtime('There are {} samples in this run. '
                  'Running off-hours module with the following parameters:\n'
                  'MiSeqPath: {},\n'
                  'MiSeqFolder: {},\n'
                  'SampleSheet: {}'.format(self.samplecount, self.miseqpath, self.miseqfolder, self.samplesheet),
                  self.start)
        # Run the fastqmover module now that the number of sequences is known
        self.fastqlinker()

    def fastqlinker(self):
        """Ensure that the sequencing run is completed, and then copy the directory to :self.path"""
        # Module-specific imports
        import time
        import re
        import shutil
        # Glob for .gz files in the appropriate subfolder of :miseqfolder. Discard 'Undetermined' files
        gzfiles = [gzfile for gzfile in glob('{}/Data/Intensities/BaseCalls/*.gz'.format(self.miseqfolder))
                   if "Undetermined" not in gzfile]
        # While loop to wait until run is complete - two .gz files are created for each sample
        while len(gzfiles) < 2 * self.samplecount:
            printtime('Waiting for run to finish. Currently, {} out of a total of {} fastq.gz files '
                      'have been created'.format(len(gzfiles), 2 * self.samplecount), self.start)
            # Sleep for five minutes
            time.sleep(300)
            # Check the number of .gz files again
            gzfiles = [gzfile for gzfile in glob('{}/Data/Intensities/BaseCalls/*.gz'.format(self.miseqfolder))
                       if "Undetermined" not in gzfile]
        # Link the .gz files to :self.path
        # Map the lambda function of symlinking each .gz file to :self.path/filename.fastq.gz
        map(lambda x: os.symlink(x, '{}{}'.format(self.path, os.path.basename(x)))
            # Don't link if the file is already present
            if not os.path.isfile('{}{}'.format(self.path, os.path.basename(x))) and
            # Don't link if a folder with a portion of the file name is present
            # e.g. 2015-SEQ-0385_S1_L001_R1_001.fastq.gz would not be linked if a folder named 2015-SEQ-0385 was present
            not os.path.isdir('{}{}'.format(self.path, re.split("_S\d+_L001", os.path.basename(x))[0]))
            # Else x (I'm not sure what this does, or why it was required) all mapped to each entry in :gzfiles
            else x, gzfiles)
        # Copy the GenerateFASTQRunStatistics.xml, RunInfo.xml, and SampleSheet.csv to self.path
        map(lambda x: shutil.copyfile('{}/{}'.format(self.miseqfolder, x), '{}{}'.format(self.path, x))
            # Don't copy if the file is already present
            if not os.path.isfile('{}{}'.format(self.path, x)) else x,
            # List of the files of interest
            ['GenerateFASTQRunStatistics.xml', 'RunInfo.xml', 'SampleSheet.csv'])

    def __init__(self, inputobject):
        """Initialise variables"""
        import sys
        self.path = inputobject.path
        self.miseqfolder = inputobject.args['f']
        self.miseqfoldername = ""
        self.customsamplesheet = inputobject.customsamplesheet
        self.start = inputobject.starttime
        self.samplecount = 0
        self.samplesheet = ""
        try:
            self.miseqpath = os.path.join(inputobject.args['m'], "")
        except AttributeError:
            print('MiSeqPath argument is required in order to use the off-hours module. Please provide this argument '
                  'and run the script again.')
            sys.exit()
        # # Assert that provided arguments are valid
        # self.assertpathsandfiles()
        # # Determine the number of samples to process by parsing the sample sheet
        # self.numberofsamples()
