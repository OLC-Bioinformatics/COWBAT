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
        self.forwardlength = self.metadata.samples.forwardlength if self.forwardlength.lower()\
            == 'full' else self.forwardlength
        # Set :reverselength to :header.reverselength
        self.reverselength = self.metadata.samples.reverselength if self.reverselength.lower() \
            == 'full' else self.reverselength
        # As the number of cycles required is the number of forward reads + the index(8) + the second index(8)
        # Also set the basemask variable as required
        if self.reverselength != '0':
            readsneeded = int(self.forwardlength) + int(self.reverselength)
            basemask = "Y{}n*,I8,I8,Y{}n*".format(self.forwardlength, self.reverselength)
            nohup = "nohup make -j 16"
        else:
            readsneeded = int(self.forwardlength) + 16
            basemask = "Y{}n*,I8,I8,n*".format(self.forwardlength)
            nohup = "nohup make -j 16 r1"
        printtime('There are {} samples in this run. '
                  'Running fastq creating module with the following parameters:\n'
                  'MiSeqPath: {},\n'
                  'MiSeqFolder: {},\n'
                  'SampleSheet: {}'.format(samplecount / 2, self.miseqpath, self.miseqfolder,
                                           '{}/SampleSheet_modified.csv'.format(self.fastqdestination)), self.start)
        # Count the number of completed cycles in the run of interest
        cycles = glob('{}Data/Intensities/BaseCalls/L001/C*'.format(self.miseqfolder))
        while len(cycles) < readsneeded:
            printtime('Currently at {} cycles. Waiting until the MiSeq reaches cycle {}'.format(len(cycles),
                      readsneeded), self.start)
            sleep(30)
            cycles = glob('{}Data/Intensities/BaseCalls/L001/C*'.format(self.miseqfolder))
        # Define the bcl2fastq system call
        bclcall = "configureBclToFastq.pl --input-dir {}Data/Intensities/BaseCalls/ " \
                  "--output-dir {} --force --sample-sheet {}/SampleSheet_modified.csv " \
                  "--mismatches 1 --no-eamss --fastq-cluster-count 0 --compression none --use-bases-mask {}"\
            .format(self.miseqfolder, self.fastqdestination, self.fastqdestination, basemask)
        # Define the nohup system call
        nohupcall = "cd {} && {}".format(self.fastqdestination, nohup)
        if not os.path.isdir("{}/Project_{}".format(self.fastqdestination, self.projectname)):
            # Call configureBclToFastq.pl
            printtime('Running bcl2fastq', self.start)
            # Run the commands
            execute(bclcall, "", 1001)
            execute(nohupcall, '{}/nohup.out'.format(self.fastqdestination), 101)
        # Populate the metadata
        for sample in self.metadata.samples:
            sample.commands.nohupcall = nohupcall
            sample.commands.bclcall = bclcall
        # Link the fastq files to a central folder so they can be processed
        self.fastqmover()

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
            print strainfastqfiles
            # Don't link the files if they have already been linked
            if strainfastqfiles < self.numreads:
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
        self.miseqout = ""
        self.projectname = 'fastqCreation'
        self.projectpath = ""
        self.numreads = inputobject.numreads
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
        self.customsamplesheet = self.assertions.customsamplesheet
        # Parse the sample sheet and other metadata files here
        self.metadata = runMetadata.Metadata(inputobject)
        self.metadata.parseruninfo()
        self.metadata.parsesamplesheet()
        # Create variables from this method
        self.flowcell = self.metadata.flowcell
        self.instrument = self.metadata.instrument
        self.samples = self.metadata.samples
        # self.header = self.metadata.
        self.ids = self.metadata.ids
        self.date = self.metadata.date
        # import json
        # print json.dumps(self.metadata., sort_keys=True, indent=4, separators=(',', ': '))
        # Create fastq files
        self.createfastq()
