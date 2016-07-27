#!/usr/bin/env python
import subprocess

import depth
import fastqCreator
import metadataprinter
import offhours
import quality
import runMetadata
import spadesRun
from accessoryFunctions import *

__author__ = 'adamkoziol'


class RunSpades(object):
    def assembly(self):
        """Helper function for file creation (if desired), manipulation, quality assessment,
        and trimming as well as the assembly"""
        import fastqmover
        # Run the fastq creation script - if argument is provided
        if self.fastqcreation:
            self.runmetadata = fastqCreator.CreateFastq(self)
            # Print the metadata to file
            metadataprinter.MetadataPrinter(self)
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
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    def quality(self):
        import errorcorrection
        import mMLST
        import quaster
        import prodigal
        import mash
        """Creates quality objects and runs the quality assessment (FastQC), and quality trimming (bbduk) on the
        supplied sequences"""
        # Run FastQC on the unprocessed fastq files
        self.qualityobject.fastqcthreader('Raw')
        # Perform quality trimming and FastQC on the trimmed files
        self.qualityobject.trimquality()
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
        # Run the error correction module
        errorcorrection.Correct(self)
        metadataprinter.MetadataPrinter(self)
        # Run FastQC on the unprocessed fastq files
        self.qualityobject.fastqcthreader('trimmedcorrected')
        # Exit if only pre-processing of data is requested
        if self.preprocess:
            printtime('Pre-processing complete', starttime)
            quit()
        # Run spades
        spadesRun.Spades(self)
        prodigal.Prodigal(self)
        metadataprinter.MetadataPrinter(self)
        # Run mash
        mash.Mash(self, 'mash')
        metadataprinter.MetadataPrinter(self)
        # Run rMLST
        mMLST.PipelineInit(self, 'rmlst')
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
        # Run quast assembly metrics
        quaster.Quast(self)
        # Calculate the depth of coverage as well as other quality metrics using Qualimap
        metadataprinter.MetadataPrinter(self)
        depth.QualiMap(self)
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    def typing(self):
        import mMLST
        import GeneSeekr
        import sixteenS
        import univec
        import prophages
        import plasmidfinder
        import serotype
        import virulence
        import armi
        import vtyper
        import coregenome
        import resfinder
        # Run modules and print metadata to file
        mMLST.PipelineInit(self, 'mlst')
        metadataprinter.MetadataPrinter(self)
        geneseekr = GeneSeekr.PipelineInit(self, 'geneseekr', True, 50)
        GeneSeekr.GeneSeekr(geneseekr)
        metadataprinter.MetadataPrinter(self)
        sixteens = GeneSeekr.PipelineInit(self, '16s', False, 95)
        sixteenS.SixteenS(sixteens)
        metadataprinter.MetadataPrinter(self)
        uni = univec.PipelineInit(self, 'univec', False, 80)
        univec.Univec(uni)
        metadataprinter.MetadataPrinter(self)
        pro = GeneSeekr.PipelineInit(self, 'prophages', False, 80)
        prophages.Prophages(pro)
        metadataprinter.MetadataPrinter(self)
        plasmid = GeneSeekr.PipelineInit(self, 'plasmidfinder', False, 80)
        plasmidfinder.PlasmidFinder(plasmid)
        metadataprinter.MetadataPrinter(self)
        sero = GeneSeekr.PipelineInit(self, 'serotype', True, 95)
        serotype.Serotype(sero)
        metadataprinter.MetadataPrinter(self)
        vir = GeneSeekr.PipelineInit(self, 'virulence', True, 70)
        virulence.Virulence(vir)
        metadataprinter.MetadataPrinter(self)
        armiobject = GeneSeekr.PipelineInit(self, 'ARMI', False, 70)
        armi.ARMI(armiobject)
        metadataprinter.MetadataPrinter(self)
        vtyper.Vtyper(self, 'vtyper')
        metadataprinter.MetadataPrinter(self)
        core = GeneSeekr.PipelineInit(self, 'coregenome', True, 70)
        coregenome.CoreGenome(core)
        metadataprinter.MetadataPrinter(self)
        res = resfinder.PipelineInit(self, 'resfinder', False, 80)
        resfinder.ResFinder(res)
        metadataprinter.MetadataPrinter(self)

    # TODO Sistr as a module

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        import reporter
        import compress
        import multiprocessing
        import versions
        printtime('Welcome to the CFIA de novo bacterial assembly pipeline {}'.format(pipelinecommit), startingtime)
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path, '')
        self.numreads = args.numreads
        self.offhours = args.offhours
        self.fastqcreation = args.fastqcreation
        self.fastqdestination = args.destinationfastq
        self.reffilepath = os.path.join(args.referencefilepath, '')
        self.forwardlength = args.readlengthforward
        self.reverselength = args.readlengthreverse
        self.numreads = 1 if self.reverselength == 0 else 2
        self.kmers = args.kmerrange
        self.preprocess = args.preprocess
        self.updatedatabases = args.updatedatabases
        # Define the start time
        self.starttime = startingtime
        self.customsamplesheet = args.customsamplesheet
        if self.customsamplesheet:
            assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}'\
                .format(self.customsamplesheet)

        self.basicassembly = args.basicassembly
        if not self.customsamplesheet and not os.path.isfile("{}SampleSheet.csv".format(self.path)):
            self.basicassembly = True
            printtime('Could not find a sample sheet. Performing basic assembly (no run metadata captured)',
                      self.starttime)
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
        # Assertions to ensure that the provided variables are valid
        make_path(self.path)
        assert os.path.isdir(self.path), u'Supplied path location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = '{}reports'.format(self.path)
        assert os.path.isdir(self.reffilepath), u'Reference file path is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = str(pipelinecommit)
        self.homepath = scriptpath
        self.runinfo = ''
        # Initialise the metadata object
        self.runmetadata = ''
        try:
            # Start the assembly
            self.assembly()
            # Create the quality object
            self.qualityobject = quality.Quality(self)
            self.quality()
            # Print the metadata to file
            metadataprinter.MetadataPrinter(self)
            # Perform typing of assemblies
            self.typing()
        except KeyboardInterrupt:
            raise KeyboardInterruptError
        # Create a report
        reporter.Reporter(self)
        compress.Compress(self)
        # Get all the versions of the software used
        versions.Versions(self)
        metadataprinter.MetadataPrinter(self)

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
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
    parser.add_argument('-o', '--offhours',
                        action='store_true',
                        help='Optionally run the off-hours module that will search for MiSeq runs in progress, wait '
                             'until the run is complete, and assemble the run')
    parser.add_argument('-F', '--fastqcreation',
                        action='store_true',
                        help='Optionally run the fastq creation module that will search for MiSeq runs in progress, '
                             'run bcl2fastq to create fastq files, and assemble the run')
    parser.add_argument('-d', '--destinationfastq',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             'Defaults to path/miseqfolder')
    parser.add_argument('-m', '--miseqpath',
                        help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', '--miseqfolder',
                        help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', '--readlengthforward',
                        default='full',
                        help='Length of forward reads to use. Can specify "full" to take the full length of forward '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-r2', '--readlengthreverse',
                        default='full',
                        help='Length of reverse reads to use. Can specify "full" to take the full length of reverse '
                             'reads specified on the SampleSheet. Defaults to "full"')
    parser.add_argument('-r', '--referencefilepath',
                        default='/spades_pipeline/SPAdesPipelineFiles',
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    parser.add_argument('-k', '--kmerrange',
                        default='21,33,55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 21,33,55,77,99,127')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet and name of sample sheet file '
                             'e.g. /home/name/folder/BackupSampleSheet.csv. Note that this sheet must still have the '
                             'same format of Illumina SampleSheet.csv files')
    parser.add_argument('-b', '--basicassembly',
                        action='store_true',
                        help='Performs a basic de novo assembly, and does not collect run metadata')
    parser.add_argument('-p', '--preprocess',
                        action='store_true',
                        help='Perform quality trimming and error correction only. Do not assemble the trimmed + '
                             'corrected reads')
    parser.add_argument('-u', '--updatedatabases',
                        action='store_true',
                        help='Optionally update (r)MLST databases')

    # Get the arguments into an object
    arguments = parser.parse_args()
    starttime = time()
    # Run the pipeline
    RunSpades(arguments, commit, starttime, homepath)
    printtime('Assembly and characterisation complete', starttime)
    quit()
