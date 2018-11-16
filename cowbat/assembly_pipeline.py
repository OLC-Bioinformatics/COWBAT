#!/usr/bin/env python3
from spadespipeline.typingclasses import GDCS, Resistance, Prophages, Serotype, Univec, Virulence
from accessoryFunctions.accessoryFunctions import MetadataObject, GenObject, make_path, SetupLogging
from spadespipeline.legacy_vtyper import Vtyper as LegacyVtyper
import accessoryFunctions.metadataprinter as metadataprinter
from sixteenS.sixteens_full import SixteenS as SixteensFull
import spadespipeline.runMetadata as runMetadata
from spadespipeline.basicAssembly import Basic
import spadespipeline.fastqmover as fastqmover
from spadespipeline.mobrecon import MobRecon
import spadespipeline.compress as compress
import spadespipeline.prodigal as prodigal
import spadespipeline.reporter as reporter
import spadespipeline.quality as quality
import spadespipeline.depth as depth
import spadespipeline.sistr as sistr
import spadespipeline.skesa as skesa
import spadespipeline.phix as phix
from MLSTsippr.mlst import GeneSippr as MLSTSippr
from metagenomefilter import automateCLARK
from genesippr.genesippr import GeneSippr
from geneseekr.blast import BLAST
import coreGenome.core as core
import MASHsippr.mash as mash
from argparse import ArgumentParser
import multiprocessing
from time import time
import logging
import os

__version__ = '0.4.4.1'
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
        metadataprinter.MetadataPrinter(inputobject=self)

    def helper(self):
        """Helper function for file creation (if desired), manipulation, quality assessment,
        and trimming as well as the assembly"""
        # Simple assembly without requiring accessory files (SampleSheet.csv, etc).
        if self.basicassembly:
            self.runmetadata = Basic(inputobject=self)
        else:
            # Populate the runmetadata object by parsing the SampleSheet.csv, GenerateFASTQRunStatistics.xml, and
            # RunInfo.xml files
            self.runinfo = os.path.join(self.path, 'RunInfo.xml')
            self.runmetadata = runMetadata.Metadata(passed=self)
            # Extract the flowcell ID and the instrument name if the RunInfo.xml file was provided
            self.runmetadata.parseruninfo()
            # Extract PhiX mapping information from the run
            phi = phix.PhiX(inputobject=self)
            phi.main()
            # Populate the lack of bclcall and nohup call into the metadata sheet
            for sample in self.runmetadata.samples:
                sample.commands = GenObject()
                sample.commands.nohupcall = 'NA'
                sample.commands.bclcall = 'NA'
            # Move/link the FASTQ files to strain-specific working directories
            fastqmover.FastqMover(inputobject=self)
        # Print the metadata to file
        metadataprinter.MetadataPrinter(inputobject=self)

    def create_quality_object(self):
        """
        Create the quality object
        """
        self.qualityobject = quality.Quality(inputobject=self)

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
        metadataprinter.MetadataPrinter(inputobject=self)
        if self.preprocess:
            logging.info('Pre-processing complete')
            quit()

    def fastq_validate(self):
        """
        Attempt to detect and fix issues with the FASTQ files
        """
        self.qualityobject.validate_fastq()
        metadataprinter.MetadataPrinter(inputobject=self)

    def fastqc_raw(self):
        """
        Run FastQC on the unprocessed FASTQ files
        """
        self.qualityobject.fastqcthreader(level='Raw')
        metadataprinter.MetadataPrinter(inputobject=self)

    def quality_trim(self):
        """
        Perform quality trimming and FastQC on the trimmed files
        """
        self.qualityobject.trimquality()
        metadataprinter.MetadataPrinter(inputobject=self)

    def fastqc_trimmed(self):
        """
        Run FastQC on the quality trimmed FASTQ files
        """
        self.qualityobject.fastqcthreader(level='Trimmed')
        metadataprinter.MetadataPrinter(inputobject=self)

    def error_correct(self):
        """
        Perform error correcting on the reads
        """
        self.qualityobject.error_correction()
        metadataprinter.MetadataPrinter(inputobject=self)

    def contamination_detection(self):
        """
        Calculate the levels of contamination in the reads
        """
        self.qualityobject.contamination_finder(report_path=self.reportpath)
        metadataprinter.MetadataPrinter(inputobject=self)

    def fastqc_trimmedcorrected(self):
        """
        Run FastQC on the processed fastq files
        """
        self.qualityobject.fastqcthreader(level='trimmedcorrected')
        metadataprinter.MetadataPrinter(inputobject=self)

    def assemble(self):
        """
        Assemble genomes and perform some basic quality analyses
        """
        # Assemble genomes
        self.assemble_genomes()
        # Calculate assembly metrics on raw assemblies
        self.quality_features(analysis='raw')
        # Calculate the depth of coverage as well as other quality metrics using Qualimap
        self.qualimap()
        # Calculate assembly metrics on polished assemblies
        self.quality_features(analysis='polished')
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
        assembly = skesa.Skesa(inputobject=self)
        assembly.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def qualimap(self):
        """
        Calculate the depth of coverage as well as other quality metrics using Qualimap
        """
        qual = depth.QualiMap(inputobject=self)
        qual.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def quality_features(self, analysis):
        """
        Extract features from assemblies such as total genome size, longest contig, and N50
        """
        features = quality.QualityFeatures(inputobject=self,
                                           analysis=analysis)
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
        g_qaml = quality.GenomeQAML(inputobject=self)
        g_qaml.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def clark(self):
        """
        Run CLARK metagenome analyses on the raw reads and assemblies if the system has adequate resources
        """
        # Run CLARK typing on the .fastq and .fasta files
        automateCLARK.PipelineInit(inputobject=self,
                                   extension='fasta',
                                   light=True
                                   )
        automateCLARK.PipelineInit(inputobject=self,
                                   extension='fastq',
                                   light=True)

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
        # Resistance finding - raw reads
        self.ressippr()
        # Resistance finding - assemblies
        self.resfinder()
        # Run MOB-suite
        self.mob_suite()
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
        mash.Mash(inputobject=self,
                  analysistype='mash')
        metadataprinter.MetadataPrinter(inputobject=self)

    def rmlst(self):
        """
        Run rMLST analyses
        """
        rmlst = MLSTSippr(args=self,
                          pipelinecommit=self.commit,
                          startingtime=self.starttime,
                          scriptpath=self.homepath,
                          analysistype='rMLST',
                          pipeline=True,
                          cutoff=1.0)
        rmlst.runner()
        metadataprinter.MetadataPrinter(inputobject=self)

    def sixteens(self):
        """
        Run the 16S analyses
        """
        SixteensFull(args=self,
                     pipelinecommit=self.commit,
                     startingtime=self.starttime,
                     scriptpath=self.homepath,
                     analysistype='sixteens_full',
                     cutoff=0.95)
        metadataprinter.MetadataPrinter(inputobject=self)

    def run_gdcs(self):
        """
        Determine the presence of genomically-dispersed conserved sequences for Escherichia, Listeria, and Salmonella
        strains
        """
        # Run the GDCS analysis
        GDCS(inputobject=self)
        metadataprinter.MetadataPrinter(inputobject=self)

    def genesippr(self):
        """
        Find genes of interest
        """
        GeneSippr(args=self,
                  pipelinecommit=self.commit,
                  startingtime=self.starttime,
                  scriptpath=self.homepath,
                  analysistype='genesippr',
                  cutoff=0.95,
                  pipeline=False,
                  revbait=False)
        metadataprinter.MetadataPrinter(inputobject=self)

    def mob_suite(self):
        """

        """
        mob = MobRecon(metadata=self.runmetadata.samples,
                       analysistype='mobrecon',
                       databasepath=self.reffilepath,
                       threads=self.cpus,
                       logfile=self.logfile,
                       reportpath=self.reportpath)
        mob.mob_recon()
        metadataprinter.MetadataPrinter(inputobject=self)

    def ressippr(self):
        """
        Resistance finding - raw reads
        """
        res = Resistance(args=self,
                         pipelinecommit=self.commit,
                         startingtime=self.starttime,
                         scriptpath=self.homepath,
                         analysistype='resfinder',
                         cutoff=0.7,
                         pipeline=False,
                         revbait=True)
        res.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def resfinder(self):
        """
        Resistance finding - assemblies
        """
        resfinder = BLAST(args=self,
                          analysistype='resfinder_assembled')
        resfinder.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def prophages(self, cutoff=90):
        """
        Prophage detection
        :param cutoff: cutoff value to be used in the analyses
        """
        prophages = Prophages(args=self,
                              analysistype='prophages',
                              cutoff=cutoff,
                              unique=True)
        prophages.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def univec(self):
        """
        Univec contamination search
        """
        univec = Univec(args=self,
                        analysistype='univec',
                        cutoff=80,
                        unique=True)
        univec.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def virulence(self):
        """
        Virulence gene detection
        """
        vir = Virulence(args=self,
                        pipelinecommit=self.commit,
                        startingtime=self.starttime,
                        scriptpath=self.homepath,
                        analysistype='virulence',
                        cutoff=0.95,
                        pipeline=False,
                        revbait=True)
        vir.reporter()
        metadataprinter.MetadataPrinter(inputobject=self)

    def typing(self):
        """
        Perform analyses that use genera-specific databases
        """
        # Run modules and print metadata to file
        # MLST
        self.mlst()
        # Serotyping
        self.serosippr()
        # Assembly-based vtyper
        self.legacy_vtyper()
        # Core genome calculation
        self.coregenome()
        # Sistr
        self.sistr()

    def mlst(self):
        """
         MLST analyses
        """
        mlst = MLSTSippr(args=self,
                         pipelinecommit=self.commit,
                         startingtime=self.starttime,
                         scriptpath=self.homepath,
                         analysistype='MLST',
                         cutoff=1.0,
                         pipeline=True)
        mlst.runner()
        metadataprinter.MetadataPrinter(inputobject=self)

    def serosippr(self):
        """
        Serotyping analyses
        """
        Serotype(args=self,
                 pipelinecommit=self.commit,
                 startingtime=self.starttime,
                 scriptpath=self.homepath,
                 analysistype='serosippr',
                 cutoff=0.90,
                 pipeline=True)
        metadataprinter.MetadataPrinter(inputobject=self)

    def legacy_vtyper(self):
        """
        Legacy vtyper - uses ePCR
        """
        legacy_vtyper = LegacyVtyper(inputobject=self,
                                     analysistype='legacy_vtyper')
        legacy_vtyper.vtyper()
        metadataprinter.MetadataPrinter(inputobject=self)

    def coregenome(self):
        """
        Core genome calculation
        """
        coregen = core.CoreGenome(args=self,
                                  analysistype='coregenome',
                                  genus_specific=True)
        coregen.seekr()
        core.AnnotatedCore(inputobject=self)
        metadataprinter.MetadataPrinter(inputobject=self)

    def sistr(self):
        """
        Sistr
        """
        sistr.Sistr(inputobject=self,
                    analysistype='sistr')
        metadataprinter.MetadataPrinter(inputobject=self)

    def __init__(self, args):
        """
        Initialises the variables required for this class
        :param args: list of arguments passed to the script
        """
        SetupLogging()
        logging.info('Welcome to the CFIA de novo bacterial assembly pipeline {}'
                     .format(__version__))
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
            logging.warning('Could not find a sample sheet. Performing basic assembly (no run metadata captured)')
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else multiprocessing.cpu_count() - 1
        # Assertions to ensure that the provided variables are valid
        make_path(self.path)
        assert os.path.isdir(self.path), 'Supplied path location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = os.path.join(self.path, 'reports')
        assert os.path.isdir(self.reffilepath), 'Reference file path is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = __version__
        self.homepath = args.homepath
        self.logfile = os.path.join(self.path, 'logfile')
        self.runinfo = str()
        self.pipeline = True
        self.qualityobject = MetadataObject()
        # Initialise the metadata object
        self.runmetadata = MetadataObject()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(__version__))
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
    # Run the pipeline
    pipeline = RunAssemble(arguments)
    pipeline.main()
    logging.info('Assembly and characterisation complete')
