#!/usr/bin/env python3
from spadespipeline.typingclasses import GDCS, ResFinder, Resistance, Prophages, Plasmids, PlasmidExtractor, Serotype, \
    Univec, Virulence
from accessoryFunctions.accessoryFunctions import MetadataObject, GenObject, printtime, make_path
from sixteenS.sixteens_full import SixteenS as SixteensFull
import spadespipeline.metadataprinter as metadataprinter
import spadespipeline.primer_finder_bbduk as vtyper
import spadespipeline.GeneSeekr as GeneSeekrMethod
import spadespipeline.runMetadata as runMetadata
from spadespipeline.basicAssembly import Basic
import spadespipeline.fastqmover as fastqmover
import spadespipeline.compress as compress
import spadespipeline.prodigal as prodigal
import spadespipeline.reporter as reporter
import spadespipeline.quality as quality
import spadespipeline.univec as univec
import spadespipeline.depth as depth
import spadespipeline.sistr as sistr
import spadespipeline.skesa as skesa
import spadespipeline.phix as phix
from MLSTsippr.mlst import GeneSippr as MLSTSippr
from metagenomefilter import automateCLARK
from genesippr.genesippr import GeneSippr
import coreGenome.core as core
import MASHsippr.mash as mash
from argparse import ArgumentParser
from psutil import virtual_memory
import multiprocessing
from time import time
import subprocess
import os

__author__ = 'adamkoziol'


class RunAssemble(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        # Start the assembly
        self.helper()
        # Create the quality object
        self.create_quality_object()
        # Run the quality analyses
        self.quality()
        # Perform assembly
        self.assemble()
        # Perform genus-agnostic typing
        self.agnostictyping()
        # Perform typing
        self.typing()
        # Create a report
        reporter.Reporter(self)
        # Compress or remove all large, temporary files created by the pipeline
        compress.Compress(self)
        metadataprinter.MetadataPrinter(self)

    def helper(self):
        """Helper function for file creation (if desired), manipulation, quality assessment,
        and trimming as well as the assembly"""
        # Simple assembly without requiring accessory files (SampleSheet.csv, etc).
        if self.basicassembly:
            self.runmetadata = Basic(self)
        else:
            # Populate the runmetadata object by parsing the SampleSheet.csv, GenerateFASTQRunStatistics.xml, and
            # RunInfo.xml files
            self.runinfo = os.path.join(self.path, 'RunInfo.xml')
            self.runmetadata = runMetadata.Metadata(self)
            # Extract the flowcell ID and the instrument name if the RunInfo.xml file was provided
            self.runmetadata.parseruninfo()
            # Extract PhiX mapping information from the run
            phi = phix.PhiX(self)
            phi.main()
            # Populate the lack of bclcall and nohup call into the metadata sheet
            for sample in self.runmetadata.samples:
                sample.commands = GenObject()
                sample.commands.nohupcall = 'NA'
                sample.commands.bclcall = 'NA'
            # Move/link the FASTQ files to strain-specific working directories
            fastqmover.FastqMover(self)
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    def create_quality_object(self):
        """
        Create the quality object
        """
        self.qualityobject = quality.Quality(self)

    def quality(self):
        """
        Creates quality objects and runs quality assessments and quality processes on the
        supplied sequences
        """
        # Validate that the FASTQ files are in the proper format, and that there are no issues e.g. different numbers
        # of forward and reverse reads, read length longer than quality score length, proper extension
        self.fastq_validate()
        # Run FastQC on the unprocessed fastq files
        self.fastqc_raw()
        # Perform quality trimming and FastQC on the trimmed files
        self.quality_trim()
        # Run FastQC on the trimmed files
        self.fastqc_trimmed()
        # Perform error correcting on the reads
        self.error_correct()
        # Detect contamination in the reads
        self.contamination_detection()
        # Run FastQC on the processed fastq files
        self.fastqc_trimmedcorrected()
        # Exit if only pre-processing of data is requested
        metadataprinter.MetadataPrinter(self)
        if self.preprocess:
            printtime('Pre-processing complete', self.starttime)
            quit()

    def fastq_validate(self):
        """
        Attempt to detect and fix issues with the FASTQ files
        """
        self.qualityobject.validate_fastq()
        metadataprinter.MetadataPrinter(self)

    def fastqc_raw(self):
        """
        Run FastQC on the unprocessed FASTQ files
        """
        self.qualityobject.fastqcthreader('Raw')
        metadataprinter.MetadataPrinter(self)

    def quality_trim(self):
        """
        Perform quality trimming and FastQC on the trimmed files
        """
        self.qualityobject.trimquality()
        metadataprinter.MetadataPrinter(self)

    def fastqc_trimmed(self):
        """
        Run FastQC on the quality trimmed FASTQ files
        """
        self.qualityobject.fastqcthreader('Trimmed')
        metadataprinter.MetadataPrinter(self)

    def error_correct(self):
        """
        Perform error correcting on the reads
        """
        self.qualityobject.error_correction()
        metadataprinter.MetadataPrinter(self)

    def contamination_detection(self):
        """
        Calculate the levels of contamination in the reads
        """
        self.qualityobject.contamination_finder()
        metadataprinter.MetadataPrinter(self)

    def fastqc_trimmedcorrected(self):
        """
        Run FastQC on the processed fastq files
        """
        self.qualityobject.fastqcthreader('trimmedcorrected')
        metadataprinter.MetadataPrinter(self)

    def assemble(self):
        """
        Assemble genomes and perform some basic quality analyses
        """
        # Assemble genomes
        self.assemble_genomes()
        # Calculate assembly metrics on raw assemblies
        self.quality_features('raw')
        # Calculate the depth of coverage as well as other quality metrics using Qualimap
        self.qualimap()
        # Calculate assembly metrics on polished assemblies
        self.quality_features('polished')
        # ORF detection
        self.prodigal()
        # Assembly quality determination
        self.genome_qaml()
        # CLARK analyses
        self.clark()

    def assemble_genomes(self):
        """
        Use skesa to assemble genomes
        """
        assembly = skesa.Skesa(self)
        assembly.main()
        metadataprinter.MetadataPrinter(self)

    def qualimap(self):
        """
        Calculate the depth of coverage as well as other quality metrics using Qualimap
        """
        qual = depth.QualiMap(self)
        qual.main()
        metadataprinter.MetadataPrinter(self)

    def quality_features(self, analysis):
        """
        Extract features from assemblies such as total genome size, longest contig, and N50
        """
        features = quality.QualityFeatures(self, analysis)
        features.main()
        metadataprinter.MetadataPrinter(self)

    def prodigal(self):
        """
        Use prodigal to detect open reading frames in the assemblies
        """
        prodigal.Prodigal(self)
        metadataprinter.MetadataPrinter(self)

    def genome_qaml(self):
        """
        Use GenomeQAML to determine the quality of the assemblies
        """
        g_qaml = quality.GenomeQAML(self)
        g_qaml.main()
        metadataprinter.MetadataPrinter(self)

    def clark(self):
        """
        Run CLARK metagenome analyses on the raw reads and assemblies if the system has adequate resources
        """
        # Determine the amount of physical memory in the system
        mem = virtual_memory()
        # If the total amount of memory is greater than 100GB (this could probably be lowered), run CLARK
        if mem.total >= 100000000000:
            # Run CLARK typing on the .fastq and .fasta files
            automateCLARK.PipelineInit(self)
            automateCLARK.PipelineInit(self, 'fastq')

        else:
            # Run CLARK typing on the .fastq and .fasta files
            automateCLARK.PipelineInit(self, light=True)
            automateCLARK.PipelineInit(self, 'fastq', light=True)
        metadataprinter.MetadataPrinter(self)

    def agnostictyping(self):
        """
        Perform typing that does not require the genus of the organism to be known
        """
        # Run mash
        self.mash()
        # Run rMLST
        self.rmlst()
        # Run the 16S analyses
        self.sixteens()
        # Calculate the presence/absence of GDCS
        self.run_gdcs()
        # Find genes of interest
        self.genesippr()
        # Plasmid finding
        self.plasmids()
        # Plasmid extracting
        self.plasmid_extractor()
        # Resistance finding - raw reads
        self.ressippr()
        # Resistance finding - assemblies
        self.resfinder()
        # Prophage detection
        self.prophages()
        # Univec contamination search
        self.univec()
        # Virulence
        self.virulence()

    def mash(self):
        """
        Run mash to determine closest refseq genome
        """
        mash.Mash(self, 'mash')
        metadataprinter.MetadataPrinter(self)

    def rmlst(self):
        """
        Run rMLST analyses
        """
        MLSTSippr(self, self.commit, self.starttime, self.homepath, 'rMLST', 1.0, True)
        metadataprinter.MetadataPrinter(self)

    def sixteens(self):
        """
        Run the 16S analyses
        """
        SixteensFull(self, self.commit, self.starttime, self.homepath, 'sixteens_full', 0.95)
        metadataprinter.MetadataPrinter(self)

    def run_gdcs(self):
        """
        Determine the presence of genomically-dispersed conserved sequences for Escherichia, Listeria, and Salmonella
        strains
        """
        # Run the GDCS analysis
        GDCS(self)
        metadataprinter.MetadataPrinter(self)

    def genesippr(self):
        """
        Find genes of interest
        """
        GeneSippr(self, self.commit, self.starttime, self.homepath, 'genesippr', 0.95, False, False)
        metadataprinter.MetadataPrinter(self)

    def plasmids(self):
        """
        Plasmid finding
        """
        Plasmids(self, self.commit, self.starttime, self.homepath, 'plasmidfinder', 0.8, False, True)
        metadataprinter.MetadataPrinter(self)

    def plasmid_extractor(self):
        """
        Extracts and types plasmid sequences
        """
        plasmids = PlasmidExtractor(self)
        plasmids.main()
        metadataprinter.MetadataPrinter(self)

    def ressippr(self):
        """
        Resistance finding - raw reads
        """
        res = Resistance(self, self.commit, self.starttime, self.homepath, 'resfinder', 0.7, False, True)
        res.main()
        metadataprinter.MetadataPrinter(self)

    def resfinder(self):
        """
        Resistance finding - assemblies
        """
        ResFinder(self)
        metadataprinter.MetadataPrinter(self)

    def prophages(self, cutoff=90):
        """
        Prophage detection
        :param cutoff: cutoff value to be used in the analyses
        """
        pro = GeneSeekrMethod.PipelineInit(self, 'prophages', False, cutoff, True)
        Prophages(pro)
        metadataprinter.MetadataPrinter(self)

    def univec(self):
        """
        Univec contamination search
        """
        uni = univec.PipelineInit(self, 'univec', False, 80, True)
        Univec(uni)
        metadataprinter.MetadataPrinter(self)

    def virulence(self):
        """
        Virulence gene detection
        """
        vir = Virulence(self, self.commit, self.starttime, self.homepath, 'virulence', 0.95, False, True)
        vir.reporter()
        metadataprinter.MetadataPrinter(self)

    def typing(self):
        """
        Perform analyses that use genera-specific databases
        """
        # Run modules and print metadata to file
        # MLST
        self.mlst()
        # Serotyping
        self.serosippr()
        # Virulence typing
        self.vtyper()
        # Core genome calculation
        self.coregenome()
        # Sistr
        self.sistr()

    def mlst(self):
        """
         MLST analyses
        """
        MLSTSippr(self, self.commit, self.starttime, self.homepath, 'MLST', 1.0, True)
        metadataprinter.MetadataPrinter(self)

    def serosippr(self):
        """
        Serotyping analyses
        """
        Serotype(self, self.commit, self.starttime, self.homepath, 'serosippr', 0.90, True)
        metadataprinter.MetadataPrinter(self)

    def vtyper(self):
        """
        Virulence typing
        """
        vtype = vtyper.PrimerFinder(self, 'vtyper')
        vtype.main()
        metadataprinter.MetadataPrinter(self)

    def coregenome(self):
        """
        Core genome calculation
        """
        coregen = GeneSeekrMethod.PipelineInit(self, 'coregenome', True, 70, False)
        core.CoreGenome(coregen)
        core.AnnotatedCore(self)
        metadataprinter.MetadataPrinter(self)

    def sistr(self):
        """
        Sistr
        """
        sistr.Sistr(self, 'sistr')
        metadataprinter.MetadataPrinter(self)

    def __init__(self, args):
        """
        Initialises the variables required for this class
        :param args: list of arguments passed to the script
        """
        printtime('Welcome to the CFIA de novo bacterial assembly pipeline {}'
                  .format(args.commit.decode('utf-8')), args.startingtime, '\033[1;94m')
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.sequencepath)
        self.reffilepath = os.path.join(args.referencefilepath)
        self.numreads = args.numreads
        self.preprocess = args.preprocess
        # Define the start time
        self.starttime = args.startingtime
        self.customsamplesheet = args.customsamplesheet
        if self.customsamplesheet:
            assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}'\
                .format(self.customsamplesheet)
        self.basicassembly = args.basicassembly
        if not self.customsamplesheet and not os.path.isfile(os.path.join(self.path, 'SampleSheet.csv')):
            self.basicassembly = True
            printtime('Could not find a sample sheet. Performing basic assembly (no run metadata captured)',
                      self.starttime)
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else multiprocessing.cpu_count() - 1
        # Assertions to ensure that the provided variables are valid
        make_path(self.path)
        assert os.path.isdir(self.path), 'Supplied path location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.reffilepath), 'Reference file path is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = args.commit.decode('utf-8')
        self.homepath = args.homepath
        self.logfile = os.path.join(self.path, 'logfile')
        self.runinfo = str()
        self.pipeline = True
        self.qualityobject = MetadataObject()
        # Initialise the metadata object
        self.runmetadata = MetadataObject()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to folder containing sequencing reads')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    parser.add_argument('-n', '--numreads',
                        default=2,
                        type=int,
                        help='Specify the number of reads. Paired-reads:'
                        ' 2, unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet and name of sample sheet file '
                             'e.g. /home/name/folder/BackupSampleSheet.csv. Note that this sheet must still have the '
                             'same format of Illumina SampleSheet.csv files')
    parser.add_argument('-b', '--basicassembly',
                        action='store_true',
                        help='Performs a basic de novo assembly, and does not collect run metadata')
    parser.add_argument('-p', '--preprocess',
                        action='store_true',
                        help='Performs quality trimming and error correction only. Does not assemble the trimmed + '
                             'corrected reads')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.startingtime = time()
    arguments.homepath = homepath
    arguments.commit = commit
    # Run the pipeline
    pipeline = RunAssemble(arguments)
    pipeline.main()
    printtime('Assembly and characterisation complete', arguments.startingtime)
