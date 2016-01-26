#!/usr/bin/env python
import os
import sys
from offhours import Offhours
from glob import glob
import runMetadata
from accessoryFunctions import make_path, printtime, execute
__author__ = 'adamkoziol'


class CreateFastq(object):

    def createfastq(self):
        """Uses bcl2fastq to create .fastq files from a MiSeqRun"""
        from time import sleep
        # Initialise samplecount
        samplecount = 0
        # If the fastq destination folder is not provided, make the default value of :path/:miseqfoldername
        self.fastqdestination = self.fastqdestination if self.fastqdestination else self.path + self.miseqfoldername
        # Make the path
        make_path(self.fastqdestination)
        # bcl2fastq requires an older version of the sample sheet, this recreates the required version
        # Create the new sample sheet
        with open('{}/SampleSheet_modified.csv'.format(self.fastqdestination), "wb") as modifiedsamplesheet:
            # Write the required headings to the file
            modifiedsamplesheet.write(
                "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n")
            for strain in self.samples:
                # Create a combined index of index1-index2
                modifiedindex = '{}-{}'.format(strain.run.index1, strain.run.index2)
                # The list of items to print to each line of the modified sample sheet
                printlist = [self.flowcell, '1', strain.name, str(strain.run.SampleNumber), modifiedindex,
                             strain.run.Description, 'N', 'NA',
                             strain.run.investigator, self.projectname]
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
            self.readsneeded = int(self.forwardlength) + int(self.reverselength) + 16
            basemask = "Y{}n*,I8,I8,Y{}n*".format(self.forwardlength, self.reverselength)
            nohup = "nohup make -j 16"
        else:
            self.readsneeded = int(self.forwardlength) + 16
            basemask = "Y{}n*,I8,I8,n*".format(self.forwardlength)
            nohup = "nohup make -j 16 r1"
        printtime('There are {} samples in this run. '
                  'Running fastq creating module with the following parameters:\n'
                  'MiSeqPath: {},\n'
                  'MiSeqFolder: {},\n'
                  'SampleSheet: {}'.format(samplecount, self.miseqpath, self.miseqfolder,
                                           '{}/SampleSheet_modified.csv'.format(self.fastqdestination)), self.start)
        # Count the number of completed cycles in the run of interest
        cycles = glob('{}Data/Intensities/BaseCalls/L001/C*'.format(self.miseqfolder))
        while len(cycles) < self.readsneeded:
            printtime('Currently at {} cycles. Waiting until the MiSeq reaches cycle {}'.format(len(cycles),
                      self.readsneeded), self.start)
            sleep(30)
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
        if not os.path.isdir("{}/Project_{}".format(self.fastqdestination, self.projectname)):
            # Call configureBclToFastq.pl
            printtime('Running bcl2fastq', self.start)
            # Run the commands
            execute(bclcall, "")
            execute(nohupcall, '{}/nohup.out'.format(self.fastqdestination))
        # Populate the metadata
        for sample in self.metadata.samples:
            sample.commands.nohupcall = nohupcall
            sample.commands.bclcall = bclcall
        # Link the fastq files to a central folder so they can be processed
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
        from re import sub
        import errno
        # Create the project path variable
        self.projectpath = self.fastqdestination + "/Project_" + self.projectname
        # Iterate through all the sample names
        for sample in self.metadata.samples:
            # Make the outputdir variable
            outputdir = '{}{}'.format(self.path, sample.name)
            # Find any fastq files with the sample name
            strainfastqfiles = glob('{}/{}*.fastq*'.format(outputdir, sample.name))
            # Don't link the files if they have already been linked
            if len(strainfastqfiles) < self.numreads:
                # Glob all the .gz files in the subfolders - projectpath/Sample_:sample.name/*.gz
                for fastq in sorted(glob('{}/Sample_{}/*.gz'.format(self.projectpath, sample.name))):
                    # Try/except loop link .gz files to self.path
                    try:
                        # Symlink fastq file to the path, but renames them first using the sample number.
                        # 2015-SEQ-1283_GGACTCCT-GCGTAAGA_L001_R1_001.fastq.gz is renamed:
                        # 2015-SEQ-1283_S1_L001_R1_001.fastq.gz
                        make_path(outputdir)
                        os.symlink(fastq, '{}/{}'.format(outputdir, os.path.basename(sub('\w{8}-\w{8}',
                                                                                         'S{}'.format(
                                                                                                 sample.run.SampleNumber
                                                                                         ), fastq))))
                    # Except os errors
                    except OSError as exception:
                        # If there is an exception other than the file exists, raise it
                        if exception.errno != errno.EEXIST:
                            raise
                # Repopulate :strainfastqfiles with the freshly-linked files
                strainfastqfiles = glob('{}/{}*.fastq*'.format(outputdir, sample.name))
            # Populate the metadata object with the name/path of the fastq files
            sample.general.fastqfiles = strainfastqfiles
            # Save the outputdir to the metadata object
            sample.run.outputdirectory = outputdir

    def __init__(self, inputobject):
        """Initialise variables"""
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.forwardlength = inputobject.forwardlength
        self.reverselength = inputobject.reverselength
        self.fastqdestination = inputobject.fastqdestination
        self.homepath = inputobject.homepath
        self.miseqout = ""
        self.projectname = 'fastqCreation'
        self.projectpath = ""
        self.numreads = inputobject.numreads
        self.readsneeded = 0
        self.commit = inputobject.commit
        try:
            self.miseqpath = os.path.join(inputobject.args['m'], "")
        except AttributeError:
            print('MiSeqPath argument is required in order to use the fastq creation module. Please provide this '
                  'argument and run the script again.')
            sys.exit()
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
