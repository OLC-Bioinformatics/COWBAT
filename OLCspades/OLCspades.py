#!/usr/bin/env python
import os
import subprocess
import runMetadata
import offhours
import fastqCreator
# import json
import metadataprinter
import spadesRun
__author__ = 'adamkoziol'


class RunSpades(object):
    def assembly(self):
        """Helper function for file creation (if desired), manipulation, quality assessment,
        and trimming as well as the assembly"""
        from accessoryFunctions import GenObject
        import fastqmover
        # Run the fastq creation script - if argument is provided
        if self.fastqcreation:
            self.runmetadata = fastqCreator.CreateFastq(self)
        # Simple assembly without requiring accessory files (SampleSheet.csv, etc).
        elif self.basicassembly:
            from basicAssembly import Basic
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

    def typing(self):
        import mMLST
        # blaster(path, cutoff, sequencepath, allelepath, organismpath, scheme, organism)
        # mMLST.blaster(self.path, 98, self, '/spades_pipeline/SPAdesPipelineFiles/rMLST', '', '', '')
        mMLST.PipelineInit(self, 'rmlst')
        # for sample in self.runmetadata.samples:
        #     print sample.mlst.datastore

    # TODO Dictreader - tsv to dictionary

    # TODO SPAdes as library
    # TODO quast as library
    # TODO Figure out what to do about GeneMark license keys
    """
    Running GeneMark...
    WARNING: License period for GeneMark has ended!
    To update license, please visit http://topaz.gatech.edu/license_download.cgi page and fill in the form.
    You should choose GeneMarkS tool and your operating system (note that GeneMark is free for non-commercial use).
    Download the license key and replace your ~/.gm_key with the updated version. After that you can restart QUAST.
    """
    """WARNING: Can't draw plots: please install python-matplotlib."""

    # TODO CGE

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path, '')
        # self.numreads = args['n']
        self.numreads = args.numreads
        self.offhours = args.offhours
        self.fastqcreation = args.fastqcreation
        self.fastqdestination = args.destinationfastq
        self.reffilepath = os.path.join(args.referencefilepath, '')
        self.forwardlength = args.readlengthforward
        self.reverselength = args.readlengthreverse
        self.numreads = 1 if self.reverselength == 0 else 2
        self.kmers = args.kmerrange
        self.customsamplesheet = args.customsamplesheet
        self.basicassembly = args.basicassembly
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else int(subprocess.Popen("awk '/^processor/ { N++} END { print N }' "
                                                          "/proc/cpuinfo", shell=True, stdout=subprocess.PIPE)
                                                          .communicate()[0].rstrip())
        # Assertions to ensure that the provided variables are valid
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = '{}reports'.format(self.path)
        assert os.path.isdir(self.reffilepath), u'Reference file path is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = pipelinecommit
        self.homepath = scriptpath
        self.runinfo = ''
        # Initialise the metadata object
        self.runmetadata = ''
        # Define the start time
        self.starttime = startingtime
        # Start the assembly
        self.assembly()
        # Run the quality trimming module
        self.quality()
        # Run spades
        spadesRun.Spades(self)
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
        # Perform typing of assemblies
        self.typing()
        # import json
        # print json.dumps([x.dump() for x in self.runmetadata.samples],
        #                  sort_keys=True, indent=4, separators=(',', ': '))


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    from accessoryFunctions import printtime
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
    parser.add_argument('-n', '--numreads', default=2, type=int, help='Specify the number of reads. Paired-reads:'
                        ' 2, unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', '--threads', help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-o', '--offhours', action='store_true', help='Optionally run the off-hours module that will '
                        'search for MiSeq runs in progress, wait until the run is complete, and assemble the run')
    parser.add_argument('-F', '--fastqcreation', action='store_true', help='Optionally run the fastq creation module'
                        'that will search for MiSeq runs in progress, run bcl2fastq to create fastq files, and '
                        'assemble the run')
    parser.add_argument('-d', '--destinationfastq', help='Optional folder path to store .fastq files created '
                        'using the fastqCreation module. Defaults to path/miseqfolder')
    parser.add_argument('-m', '--miseqpath', help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder', help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', '--readlengthforward', default='full', help='Length of forward reads to use. Can '
                        'specify "full" to take the full length of forward reads specified on the SampleSheet. '
                        'Defaults to full')
    parser.add_argument('-r2', '--readlengthreverse', default='full', help='Length of reverse reads to use. '
                        'Can specify "full" to take the full length of reverse reads specified on the SampleSheet. '
                        'Defaults to full')
    parser.add_argument('-r', '--referencefilepath', default="/spades_pipeline/SPAdesPipelineFiles",
                        help='Provide the location of the folder containing the pipeline accessory files '
                        '(reference genomes, MLST data, etc.')
    parser.add_argument('-k', '--kmerrange', default='21,33,55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 21,33,55,77,99,127')
    parser.add_argument('-c', '--customsamplesheet', help='Path of folder containing a custom sample '
                        'sheet and name of sample sheet file e.g. /home/name/folder/BackupSampleSheet.csv. Note that '
                        'this sheet must still have the same format of Illumina SampleSheet.csv files')
    parser.add_argument('-b', '--basicassembly', action='store_true', help='Performs a basic de novo assembly, '
                        'and does not collect metadata')

    # Get the arguments into a list
    # arguments = vars(parser.parse_args())
    arguments = parser.parse_args()
    starttime = time()
    # Run the pipeline
    RunSpades(arguments, commit, starttime, homepath)
    printtime('Assembly and characterisation complete', starttime)
