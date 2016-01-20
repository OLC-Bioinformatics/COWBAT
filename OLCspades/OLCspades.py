#!/usr/bin/env python
import os
import subprocess
import runMetadata
import offhours
import fastqCreator
# import json
import metadataprinter
__author__ = 'adamkoziol'


class KeyboardInterruptError(Exception):
    pass


class RunSpades(object):
    def filemanipulation(self):
        from accessoryFunctions import GenObject
        import fastqmover
        """Helper function for file creation (if desired), manipulation, quality assessment,
        and trimming as well as the assembly"""
        # Run the fastq creation script - if argument is provided
        if self.fastqcreation:
            self.runmetadata = fastqCreator.CreateFastq(self)
        # Simple assembly without requiring accessory files (SampleSheet.csv, etc).
        elif self.basicassembly:
            self.runmetadata = Basic(self)
        else:
            # Populate the runmetadata object by parsing the SampleSheet.csv, GenerateFASTQRunStatistics.xml, and
            # RunInfo.xml files
            self.runinfo = "{}RunInfo.xml".format(self.path)
            self.runmetadata = runMetadata.Metadata(self)
            # Extract the flowcell ID and the instrument name if the RunInfo.xml file was provided
            self.runmetadata.parseruninfo()
            # Populate the lack of bclcall and nohup call into the metadata sheet
            for sample in self.runmetadata.samples:
                sample.commands = GenObject()
                sample.commands.nohupcall = 'NA'
                sample.commands.bclcall = 'NA'
            # Run the offhours fastq linking script - if argument
            if self.offhours:
                offhoursobj = offhours.Offhours(self)
                offhoursobj.assertpathsandfiles()
                offhoursobj.numberofsamples()
            # Move the files
            else:
                fastqmover.FastqMover(self)
        # Run the quality trimming module
        self.quality()

        # import json
        # print json.dumps([x.dump() for x in self.runmetadata.samples],
        #                  sort_keys=True, indent=4, separators=(',', ': '))

    def quality(self):
        """Creates quality objects and runs the quality assessment (FastQC), and quality trimming (bbduk) on the
        supplied sequences"""
        import quality
        # Create the quality object
        qualityobject = quality.Quality(self)
        # Run FastQC on the unprocessed fastq files
        qualityobject.fastqcthreader('Raw')
        # Perform quality trimming and FastQC on the trimmed files
        qualityobject.trimquality()
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    # TODO Dictreader - tsv to dictionary

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args['path'], "")
        self.numreads = args['n']
        self.offhours = args['offHours']
        self.fastqcreation = args['FastqCreation']
        self.fastqdestination = args['d']
        self.reffilepath = os.path.join(args['r'], "")
        self.forwardlength = args['r1']
        self.reverselength = args['r2']
        self.numreads = 1 if self.reverselength == 0 else 2
        self.kmers = args['k']
        self.customsamplesheet = args['c']
        self.basicassembly = args['basicAssembly']
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args['t'] if args['t'] else int(subprocess.Popen("awk '/^processor/ { N++} END { print N }' "
                                                                     "/proc/cpuinfo", shell=True,
                                                                     stdout=subprocess.PIPE)
                                                    .communicate()[0].rstrip())
        # Assertions to ensure that the provided variables are valid
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        assert os.path.isdir(self.reffilepath), u'Output location is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = pipelinecommit
        self.homepath = scriptpath
        self.runinfo = ""
        # Initialise the metadata object
        self.runmetadata = ""
        # Define the start time
        self.starttime = startingtime


class Basic(object):
    def basic(self):
        from accessoryFunctions import GenObject, MetadataObject, filer, make_path
        from glob import glob
        import os
        import errno
        # Grab any .fastq files in the path
        fastqfiles = glob('{}*fastq*'.format(self.path))
        # Extract the base name of the globbed name + path provided
        fastqnames = map(lambda x: os.path.split(x)[1], filer(fastqfiles))
        # Iterate through the names of the fastq files
        for fastqname in fastqnames:
            # Set the name
            metadata = MetadataObject()
            metadata.name = fastqname
            # Set the destination folder
            outputdir = '{}{}'.format(self.path, fastqname)
            # Make the destination folder
            make_path(outputdir)
            # Get the fastq files specific to the fastqname
            specificfastq = glob('{}{}*fastq*'.format(self.path, fastqname))
            # Link the files to the output folder
            try:
                # Link the .gz files to :self.path/:filename
                map(lambda x: os.symlink(x, '{}/{}'.format(outputdir, os.path.split(x)[1])), specificfastq)
            # Except os errors
            except OSError as exception:
                # If there is an exception other than the file exists, raise it
                if exception.errno != errno.EEXIST:
                    raise
            # Initialise the general category
            metadata.general = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = glob('{}/{}*fastq*'.format(outputdir, fastqname))
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir
            # Append the metadata to the list of samples
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.samples = []
        self.path = inputobject.path
        self.basic()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',  help='Specify path')
    parser.add_argument('-n', metavar='numreads', default=2, type=int, help='Specify the number of reads. Paired-reads:'
                        ' 2, unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', metavar='threads', help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-o', '--offHours', action='store_true', help='Optionally run the off-hours module that will '
                        'search for MiSeq runs in progress, wait until the run is complete, and assemble the run')
    parser.add_argument('-F', '--FastqCreation', action='store_true', help='Optionally run the fastq creation module'
                        'that will search for MiSeq runs in progress, run bcl2fastq to create fastq files, and '
                        'assemble the run')
    parser.add_argument('-d', metavar='destinationfastq', help='Optional folder path to store .fastq files created '
                        'using the fastqCreation module. Defaults to path/miseqfolder')
    parser.add_argument('-m', metavar='miSeqPath', help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', metavar='miseqfolder', help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', metavar='readLengthForward', default='full', help='Length of forward reads to use. Can '
                        'specify "full" to take the full length of forward reads specified on the SampleSheet. '
                        'Defaults to full')
    parser.add_argument('-r2', metavar='readLengthReverse', default='full', help='Length of reverse reads to use. '
                        'Can specify "full" to take the full length of reverse reads specified on the SampleSheet. '
                        'Defaults to full')
    # parser.add_argument('-P', metavar='projectName', help='A name for the analyses. If nothing is provided, then '
    #                     'the "Sample_Project" field in the provided sample sheet will be used. Please note that '
    #                     'bcl2fastq creates subfolders using the project name, so if multiple names are provided, the '
    #                     'results will be split as into multiple projects')
    parser.add_argument('-r', metavar='referenceFilePath', default="/spades_pipeline/SPAdesPipelineFiles",
                        help='Provide the location of the folder containing the pipeline accessory files '
                        '(reference genomes, MLST data, etc.')
    parser.add_argument('-k', metavar='kmerRange', default='21,33,55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 21,33,55,77,99,127')
    parser.add_argument('-c', metavar='customSampleSheet', help='Path of folder containing a custom sample '
                        'sheet and name of sample sheet file e.g. /home/name/folder/BackupSampleSheet.csv. Note that '
                        'this sheet must still have the same format of Illumina SampleSheet.csv files')
    parser.add_argument('-b', '--basicAssembly', action='store_true', help='Performs a basic de novo assembly, '
                        'and does not collect metadata')

    # Get the arguments into a list
    arguments = vars(parser.parse_args())

    starttime = time()
    # Run the pipeline
    output = RunSpades(arguments, commit, starttime, homepath)
    output.filemanipulation()
    print "\nElapsed Time: {:.2f} seconds".format(time() - starttime)
