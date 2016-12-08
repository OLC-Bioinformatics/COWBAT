#!/usr/bin/env python
from glob import glob
import runMetadata
from accessoryFunctions import *
from offhours import Offhours

__author__ = 'adamkoziol'


class CreateFastq(object):

    def createfastq(self):
        """Uses bcl2fastq to create .fastq files from a MiSeqRun"""
        from time import sleep
        from subprocess import call
        # Initialise samplecount
        samplecount = 0
        # If the fastq destination folder is not provided, make the default value of :path/:miseqfoldername
        self.fastqdestination = self.fastqdestination if self.fastqdestination else self.path + self.miseqfoldername
        # Make the path
        make_path(self.fastqdestination)
        # Initialise variables for storing index information
        index = ''
        indexlength = int()
        # bcl2fastq requires an older version of the sample sheet, this recreates the required version
        # Create the new sample sheet
        with open('{}/SampleSheet_modified.csv'.format(self.fastqdestination), "wb") as modifiedsamplesheet:
            # Write the required headings to the file
            modifiedsamplesheet.write(
                "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n")
            for strain in self.samples:
                # Create a combined index of index1-index2
                try:
                    strain.run.modifiedindex = '{}-{}'.format(strain.run.index, strain.run.index2)
                    indexlength = 16
                    index = 'I8,I8'
                except KeyError:
                    strain.run.modifiedindex = strain.run.index
                    indexlength = 6
                    index = 'I6'
                # The list of items to print to each line of the modified sample sheet
                printlist = [self.flowcell, '1', strain.name, str(strain.run.SampleNumber), strain.run.modifiedindex,
                             strain.run.Description, 'N', 'NA',
                             strain.run.InvestigatorName, self.projectname]
                modifiedsamplesheet.write('{}\n'.format(",".join(printlist)))
                samplecount += 1
        # Set :forward/reverse length to :header.forward/reverse length if the argument is not provided, or it's 'full',
        # otherwise  use the supplied argument
        self.forwardlength = self.metadata.header.forwardlength if self.forwardlength.lower()\
            == 'full' else self.forwardlength
        # Set :reverselength to :header.reverselength
        self.reverselength = self.metadata.header.reverselength if self.reverselength.lower() \
            == 'full' else self.reverselength
        # As the number of cycles required is the number of forward reads + the index(8) + the second index(8)
        # Also set the basemask variable as required
        if self.reverselength != '0':
            self.readsneeded = int(self.forwardlength) + int(self.reverselength) + indexlength
            basemask = "Y{}n*,{},Y{}n*".format(self.forwardlength, index, self.reverselength)
            nohup = "nohup make -j 16 > nohup.out"
        else:
            #  + 1
            self.readsneeded = int(self.forwardlength) + indexlength
            basemask = "Y{}n*,{},n*".format(self.forwardlength, index)
            nohup = "nohup make -j 16 r1 > nohup.out"
        # Handle plurality appropriately
        samples = 'samples' if samplecount > 1 else 'sample'
        number = 'are' if samplecount > 1 else 'is'
        printtime('There {} {} {} in this run. '
                  'Running fastq creating module with the following parameters:\n'
                  'MiSeqPath: {},\n'
                  'MiSeqFolder: {},\n'
                  'Fastq destination: {},\n'
                  'SampleSheet: {}'
                  .format(number, samplecount, samples, self.miseqpath, self.miseqfolder,
                          self.fastqdestination, '{}/SampleSheet_modified.csv'.format(self.fastqdestination)),
                  self.start)
        # Count the number of completed cycles in the run of interest
        cycles = glob('{}Data/Intensities/BaseCalls/L001/C*'.format(self.miseqfolder))
        while len(cycles) < self.readsneeded:
            printtime('Currently at {} cycles. Waiting until the MiSeq reaches cycle {}'.format(len(cycles),
                      self.readsneeded), self.start)
            sleep(1800)
            cycles = glob('{}Data/Intensities/BaseCalls/L001/C*'.format(self.miseqfolder))
        # configureBClToFastq requires :self.miseqfolder//Data/Intensities/BaseCalls/config.xml in order to work
        # When you download runs from BaseSpace, this file is not provided. There is an empty config.xml file that
        # can be populated with run-specific values and moved to the appropriate folder
        if not os.path.isfile('{}Data/Intensities/BaseCalls/config.xml'.format(self.miseqfolder)):
            self.configfilepopulator()
        # Define the bcl2fastq system call
        bclcall = "configureBclToFastq.pl --input-dir {}Data/Intensities/BaseCalls " \
                  "--output-dir {} --force --sample-sheet {}/SampleSheet_modified.csv " \
                  "--mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask {}"\
            .format(self.miseqfolder, self.fastqdestination, self.fastqdestination, basemask)
        # Define the nohup system call
        nohupcall = "cd {} && {}".format(self.fastqdestination, nohup)
        fnull = open(os.devnull, 'wb')
        if not os.path.isdir("{}/Project_{}".format(self.fastqdestination, self.projectname)):
            # Call configureBclToFastq.pl
            printtime('Running bcl2fastq', self.start)
            # Run the commands
            call(bclcall, shell=True, stdout=fnull, stderr=fnull)
            call(nohupcall, shell=True, stdout=fnull, stderr=fnull)
        # Populate the metadata
        for sample in self.metadata.samples:
            sample.commands = GenObject()
            sample.commands.nohup = nohupcall
            sample.commands.bcl = bclcall
            sample.run.forwardlength = self.forwardlength
            sample.run.reverselength = self.reverselength
        # Copy the fastq files to a central folder so they can be processed
        self.fastqmover()

    def configfilepopulator(self):
        """Populates an unpopulated config.xml file with run-specific values and creates
        the file in the appropriate location"""
        # Import ElementTree - try first to import the faster C version, if that doesn't
        # work, try to import the regular version
        try:
            import xml.etree.cElementTree as ElementTree
        except ImportError:
            import xml.etree.ElementTree as ElementTree
        # Set the number of cycles for each read and index using the number of reads specified in the sample sheet
        self.forwardlength = self.metadata.header.forwardlength
        self.reverselength = self.metadata.header.reverselength
        # Create a list of lists containing [cycle start, cycle end, and :runid] for each of forward reads, index 1
        # index 2, and reverse reads
        cycles = [[1, self.forwardlength, self.runid],
                  [self.forwardlength + 1, self.forwardlength + 8, self.runid],
                  [self.forwardlength + 9, self.forwardlength + 16, self.runid],
                  [self.forwardlength + 17, self.forwardlength + 16 + self.reverselength, self.runid]]
        # A dictionary of parameters (keys) and the values to use when repopulating the config file
        parameters = {'RunFolder': self.runid, 'RunFolderDate': self.metadata.date.replace("-", ""),
                      'RunFolderId': self.metadata.runnumber, 'RunFlowcellId': self.metadata.flowcell}
        # Load the xml file using element tree
        config = ElementTree.parse("{}/config.xml".format(self.homepath))
        # Get the root of the tree
        configroot = config.getroot()
        # The run node is the only child node of the root
        for run in configroot:
            # Iterate through the child nodes. There are three nodes sections that must be populated
            for child in run:
                # Find the cycles tag
                if child.tag == 'Cycles':
                    # Set the attributes with a dictionary containing the total reads
                    child.attrib = {'Last': '{}'.format(self.forwardlength + 16 + self.reverselength),
                                    'Number': '{}'.format(self.totalreads), 'First': '1'}
                elif child.tag == 'RunParameters':
                    # Name the child as runparameter for easier coding
                    runparameters = child
                    for runparameter in runparameters:
                        # This replaces data in both 'ImagingReads' and 'Reads' nodes
                        if 'Reads' in runparameter.tag:
                            # Enumerate through the run parameters
                            for indexcount, reads in enumerate(runparameter):
                                # The values for the index are 1, 2, 3, 4. Subtract one to get the index of the first
                                # list in cycles
                                index = int(runparameter.attrib['Index']) - 1
                                # Set the text value as the appropriate value from cycles
                                reads.text = str(cycles[index][indexcount])
                        # Populate the instrument value
                        if runparameter.tag == 'Instrument':
                            runparameter.text = self.instrument
                        # Iterate through the parameters in the parameter dictionary
                        for parameter in parameters:
                            # If the key is encountered
                            if runparameter.tag == parameter:
                                # Replace the text with the value
                                runparameter.text = parameters[parameter]
                        if 'Barcode' in runparameter.tag:
                            for cycle, barcode in enumerate(runparameter):
                                # Add the barcode cycles. These are the number of forward reads (+ 1 as the barcode
                                # starts 1 cycle after the first run) plus the current iterator
                                barcode.text = str(self.forwardlength + 1 + cycle)
        # Write the modified config file to the desired location
        config.write('{}Data/Intensities/BaseCalls/config.xml'.format(self.miseqfolder))

    def fastqmover(self):
        """Links .fastq files created above to :self.path/:sample.name/"""
        from shutil import move, copyfile
        import errno
        from re import sub
        # Create the project path variable
        self.projectpath = self.fastqdestination + "/Project_" + self.projectname
        # Iterate through all the sample names
        for sample in self.metadata.samples:
            # Glob all the .gz files in the subfolders - projectpath/Sample_:sample.name/*.gz
            for fastq in sorted(glob('{}/Sample_{}/*.gz'.format(self.projectpath, sample.name))):
                # Try/except loop link .gz files to self.path
                try:
                    # Move fastq file to the path, but renames them first using the sample number.
                    move(
                        fastq, '{}{}'.format(self.path, os.path.basename(
                            sub('\w{8}-\w{8}', 'S{}'.format(
                                sample.run.SampleNumber), fastq))))
                # Except os errors
                except OSError as exception:
                    # If there is an exception other than the file exists, raise it
                    if exception.errno != errno.EEXIST:
                        raise
            # Repopulate .strainfastqfiles with the freshly-linked files
            fastqfiles = glob('{}/{}*.fastq*'.format(self.fastqdestination, sample.name))
            fastqfiles = [fastq for fastq in fastqfiles if 'trimmed' not in fastq]
            # Populate the metadata object with the name/path of the fastq files
            sample.general.fastqfiles = fastqfiles
            # Save the outputdir to the metadata object
            sample.run.outputdirectory = self.fastqdestination
        # Copy the sample sheet and the run info files to the path
        copyfile(self.assertions.samplesheet, os.path.join(self.path, 'SampleSheet.csv'))
        copyfile(os.path.join(self.miseqfolder, 'RunInfo.xml'), os.path.join(self.path, 'RunInfo.xml'))

    def __init__(self, inputobject):
        """Initialise variables"""
        self.path = os.path.join(inputobject.path, '')
        self.start = inputobject.starttime
        self.fastqdestination = inputobject.fastqdestination
        self.homepath = inputobject.homepath
        self.miseqout = ""
        self.projectname = 'fastqCreation'
        self.projectpath = ""
        self.numreads = inputobject.numreads
        self.forwardlength = inputobject.forwardlength
        self.reverselength = inputobject.reverselength if self.numreads > 1 else '0'
        self.readsneeded = 0
        self.commit = inputobject.commit
        if inputobject.miseqpath:
            self.miseqpath = os.path.join(inputobject.miseqpath, "")
        else:
            print inputobject.miseqfolder
            print('MiSeqPath argument is required in order to use the fastq creation module. Please provide this '
                  'argument and run the script again.')
            quit()
        self.customsamplesheet = inputobject.customsamplesheet
        if self.customsamplesheet:
            assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}' \
                .format(self.customsamplesheet)
        # Use the assertions module from offhours to validate whether provided arguments are valid
        self.assertions = Offhours(inputobject)
        self.assertions.assertpathsandfiles()
        # Populate variables from this object
        self.miseqfolder = self.assertions.miseqfolder
        self.miseqfoldername = self.assertions.miseqfoldername
        self.customsamplesheet = self.assertions.customsamplesheet if self.assertions.customsamplesheet \
            else '{}SampleSheet.csv'.format(self.miseqfolder)
        self.runinfo = '{}RunInfo.xml'.format(self.miseqfolder)
        # Parse the sample sheet and other metadata files here
        self.metadata = runMetadata.Metadata(self)
        self.metadata.parseruninfo()
        # Create variables from this method
        self.flowcell = self.metadata.flowcell
        self.instrument = self.metadata.instrument
        self.samples = self.metadata.samples
        self.runid = self.metadata.runid
        # self.header = self.metadata.
        self.ids = self.metadata.ids
        self.date = self.metadata.date
        self.totalreads = self.metadata.totalreads
        # Create fastq files
        self.createfastq()

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    import subprocess
    from time import time
    from accessoryFunctions import printtime
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify path')
    parser.add_argument('-n', '--numreads',
                        default=2,
                        type=int,
                        help='Specify the number of reads. Paired-reads:'
                        ' 2, unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-d', '--fastqdestination',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-m', '--miseqpath',
                        required=True,
                        help='Path of the folder containing MiSeq run data folder e.g. /mnt/MiSeq')
    parser.add_argument('-f', '--miseqfolder',
                        required=True,
                        help='Name of the folder containing MiSeq run data e.g. 161129_M02466_0007_000000000-AW5L5')
    parser.add_argument('-r1', '--forwardlength',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of forward '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-r2', '--reverselength',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of reverse '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet and name of sample sheet file '
                             'e.g. /home/name/folder/BackupSampleSheet.csv. Note that this sheet must still have the '
                             'same format of Illumina SampleSheet.csv files')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.starttime = time()
    arguments.commit = commit
    arguments.homepath = homepath
    # Run the pipeline
    CreateFastq(arguments)
    # Print a bold, green exit statement
    print '\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time() - arguments.starttime) + '\033[0m'
