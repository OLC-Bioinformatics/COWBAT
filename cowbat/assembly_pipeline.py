#!/usr/bin/env python3

"""
COWBAT pipeline for assembly and typing of bacterial genomes
"""

# Standard imports
import logging
import multiprocessing
import os
from argparse import ArgumentParser
from time import time

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    GenObject,
    make_path,
    MetadataObject,
    SetupLogging,
)
from olctools.accessoryFunctions import metadataprinter

from genemethods.assemblypipeline import (
    assembly_evaluation,
    basicAssembly,
    compress,
    ec_typer,
    fastqmover,
    mobrecon,
    phix,
    prodigal,
    quality,
    reporter,
    runMetadata,
    seqsero,
    sistr,
    skesa,
)
from genemethods.assemblypipeline.primer_finder_ipcress import (
    VtyperIP as IdentityVtyper,
)
from genemethods.geneseekr.blast import BLAST
from genemethods.genesippr.genesippr import GeneSippr
from genemethods.MASHsippr import mash
from genemethods.MLST.mlst_kma import KMAMLST
from genemethods.MLSTsippr.mlst import ReportParse
from genemethods.sixteenS.sixteens_full import SixteenS as SixteensFull
from genemethods.typingclasses.typingclasses import (
    GDCS,
    Prophages,
    Resistance,
    Serotype,
    Univec,
    Verotoxin,
    Virulence,
)

# Local imports
from cowbat.version import __version__
__author__ = 'adamkoziol'


class RunAssemble:
    """
    Collection of methods to run the COWBAT pipeline
    """

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
        self.agnostic_typing()
        # Perform typing
        self.typing()
        # Compress or remove all large, temporary files created by the pipeline
        if not self.debug:
            compress.Compress(self)
        metadataprinter.MetadataPrinter(inputobject=self)

    def helper(self):
        """
        Helper function for file creation (if desired), manipulation, quality
        assessment, and trimming as well as the assembly
        """
        # Simple assembly without requiring accessory files (SampleSheet.csv,
        # etc).
        if self.basic_assembly:
            self.runmetadata = basicAssembly.Basic(inputobject=self)
        else:
            # Populate the runmetadata object by parsing the SampleSheet.csv,
            # GenerateFASTQRunStatistics.xml, and RunInfo.xml files
            self.run_info = os.path.join(self.path, 'RunInfo.xml')
            self.runmetadata = runMetadata.Metadata(passed=self)
            # Extract the flowcell ID and the instrument name if the
            # RunInfo.xml file was provided
            self.runmetadata.parseruninfo()
            # Extract PhiX mapping information from the run
            phi = phix.PhiX(inputobject=self)
            phi.main()
            # Populate the lack of bclcall and nohup call into the metadata
            # sheet
            for sample in self.runmetadata.samples:
                sample.commands = GenObject()
                sample.commands.nohupcall = 'NA'
                sample.commands.bclcall = 'NA'
            # Move/link the FASTQ files to strain-specific working directories
            fastqmover.FastqMover(
                metadata=self.runmetadata,
                path=self.sequence_path
            )
        # Print the metadata to file
        metadataprinter.MetadataPrinter(inputobject=self)

    def create_quality_object(self):
        """
        Create the quality object
        """
        self.qualityobject = quality.Quality(inputobject=self)

    def quality(self):
        """
        Creates quality objects and runs quality assessments and quality
        processes on the supplied sequences
        """
        # Validate that the FASTQ files are in the proper format, and that
        # there are no issues e.g. different numbers of forward and reverse
        # reads, read length longer than quality score length, proper extension
        if not self.debug:
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
        self.fastqc_trimmed_corrected()
        # Fix issue with bbmap gzip
        self.fix_gzip()
        # Exit if only pre-processing of data is requested
        metadataprinter.MetadataPrinter(inputobject=self)
        if self.preprocess:
            logging.info('Pre-processing complete')
            raise SystemExit

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
        self.qualityobject.contamination_finder(
            report_path=self.report_path,
            debug=self.debug,
            threads=self.cpus
        )
        metadataprinter.MetadataPrinter(inputobject=self)

    def fastqc_trimmed_corrected(self):
        """
        Run FastQC on the processed fastq files
        """
        self.qualityobject.fastqcthreader(level='trimmedcorrected')
        metadataprinter.MetadataPrinter(inputobject=self)

    def fix_gzip(self):
        """
        Fix issue with BBDuk gzip
        """
        self.qualityobject.fix_gzip()
        metadataprinter.MetadataPrinter(inputobject=self)

    def assemble(self):
        """
        Assemble genomes and perform some basic quality analyses
        """
        # Assemble genomes
        self.assemble_genomes()
        # Calculate assembly metrics on raw assemblies
        self.evaluate_assemblies()
        # ORF detection
        self.prodigal()

    def assemble_genomes(self):
        """
        Use skesa to assemble genomes
        """
        assembly = skesa.Skesa(
            log_file=self.logfile,
            metadata=self.runmetadata.samples,
            report_path=self.report_path,
            sequence_path=self.path,
            threads=self.cpus
        )
        self.runmetadata.samples = assembly.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def evaluate_assemblies(self):
        """
        Evaluate assemblies with Quast
        """
        qual = assembly_evaluation.AssemblyEvaluation(
            log_file=self.logfile,
            metadata=self.runmetadata.samples,
            sequence_path=self.path,
            threads=self.cpus
        )
        self.runmetadata.samples = qual.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def prodigal(self):
        """
        Use prodigal to detect open reading frames in the assemblies
        """
        prod = prodigal.Prodigal(
            log_file=self.logfile,
            metadata=self.runmetadata.samples
        )
        
        self.runmetadata.samples = prod.main()
        metadataprinter.MetadataPrinter(self)

    def agnostic_typing(self):
        """
        Perform typing that does not require the genus of the organism to be
        known
        """
        # Run mash
        self.mash()
        # Run rMLST on assemblies
        self.rmlst_assembled()
        # Create reports summarising the run and sample qualities
        self.quality_report()
        # Run the 16S analyses
        self.sixteens()
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
        # cgMLST
        # self.cgmlst()

    def mash(self):
        """
        Run mash to determine closest refseq genome
        """
        mash.Mash(
            inputobject=self,
            analysistype='mash'
        )
        metadataprinter.MetadataPrinter(inputobject=self)

    def rmlst_assembled(self):
        """
        Run rMLST analyses on assemblies
        """
        if not os.path.isfile(os.path.join(self.report_path, 'rmlst.csv')):
            rmlst = BLAST(
                args=self,
                analysistype='rmlst',
                cutoff=100
            )
            rmlst.seekr()
        else:
            parse = ReportParse(
                args=self,
                analysistype='rmlst'
            )
            parse.report_parse()
        metadataprinter.MetadataPrinter(inputobject=self)

    def quality_report(self):
        """
        Create reports summarising the run and sample quality outputs
        """
        qual_report = reporter.Reporter(self)
        qual_report.run_quality_reporter()
        qual_report.sample_quality_report()

    def sixteens(self):
        """
        Run the 16S analyses
        """
        SixteensFull(
            args=self,
            pipelinecommit=self.commit,
            startingtime=self.start_time,
            scriptpath=self.home_path,
            analysistype='sixteens_full',
            cutoff=0.95
        )
        metadataprinter.MetadataPrinter(inputobject=self)

    def genesippr(self):
        """
        Find genes of interest
        """
        GeneSippr(
            args=self,
            pipelinecommit=self.commit,
            startingtime=self.start_time,
            scriptpath=self.home_path,
            analysistype='genesippr',
            cutoff=0.95,
            kmer_size=13,
            pipeline=False,
            revbait=False
        )
        metadataprinter.MetadataPrinter(inputobject=self)

    def mob_suite(self):
        """
        Run MOB Suite analyses
        """
        mob = mobrecon.MobRecon(
            metadata=self.runmetadata.samples,
            analysistype='mobrecon',
            databasepath=self.ref_file_path,
            threads=self.cpus,
            logfile=self.logfile,
            reportpath=self.report_path
        )
        mob.mob_recon()
        metadataprinter.MetadataPrinter(inputobject=self)

    def ressippr(self):
        """
        Resistance finding - raw reads
        """
        res = Resistance(
            args=self,
            pipelinecommit=self.commit,
            startingtime=self.start_time,
            scriptpath=self.home_path,
            analysistype='resfinder',
            cutoff=0.7,
            pipeline=False,
            revbait=True
        )
        res.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def resfinder(self):
        """
        Resistance finding - assemblies
        """
        resfinder = BLAST(
            args=self,
            analysistype='resfinder_assembled'
        )
        resfinder.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def prophages(self, cutoff=90):
        """
        Prophage detection
        :param cutoff: cutoff value to be used in the analyses
        """
        prophages = Prophages(
            args=self,
            analysistype='prophages',
            cutoff=cutoff,
            unique=True
        )
        if not os.path.isfile(os.path.join(self.report_path, 'prophages.csv')):
            prophages.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def univec(self):
        """
        Univec contamination search
        """
        if not os.path.isfile(os.path.join(self.report_path, 'univec.csv')):
            univec = Univec(
                args=self,
                analysistype='univec',
                cutoff=80,
                unique=True
            )
            univec.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def virulence(self):
        """
        Virulence gene detection
        """
        vir = Virulence(
            args=self,
            pipelinecommit=self.commit,
            startingtime=self.start_time,
            scriptpath=self.home_path,
            analysistype='virulence',
            cutoff=0.9,
            pipeline=False,
            revbait=True
        )
        if not os.path.isfile(os.path.join(self.report_path, 'virulence.csv')):
            vir.reporter()
        metadataprinter.MetadataPrinter(inputobject=self)

    def cgmlst(self):
        """
        Run rMLST analyses on raw reads
        """
        if not os.path.isfile(os.path.join(self.report_path, 'cgmlst.csv')):
            cgmlst = KMAMLST(
                args=self,
                pipeline=True,
                analysistype='cgmlst',
                cutoff=98,
                kma_kwargs=' -cge -and'
            )
            cgmlst.main()
        else:
            parse = ReportParse(args=self,
                                analysistype='cgmlst')
            parse.report_parse()
        metadataprinter.MetadataPrinter(inputobject=self)

    def typing(self):
        """
        Perform analyses that use genera-specific databases
        """
        # Run modules and print metadata to file
        # MLST on assemblies
        self.mlst_assembled()
        # Assembly-based serotyping
        self.ec_typer()
        # Serotyping
        self.serosippr()
        # SeqSero
        self.seqsero()
        # Assembly-based vtyper
        self.identity_vtyper()
        # Raw read verotoxin typing
        self.verotoxin()
        # Sistr
        self.sistr()
        # Calculate the presence/absence of GDCS
        self.run_gdcs()
        # Create a final summary report
        self.run_report()

    def mlst_assembled(self):
        """
        Run rMLST analyses on assemblies
        """
        if not os.path.isfile(os.path.join(self.report_path, 'mlst.csv')):

            mlst = BLAST(
                args=self,
                analysistype='mlst',
                cutoff=100,
                genus_specific=True
            )
            mlst.seekr()
        else:
            parse = ReportParse(args=self,
                                analysistype='mlst')
            parse.report_parse()
        metadataprinter.MetadataPrinter(inputobject=self)

    def ec_typer(self):
        """
        Assembly-based serotyping
        """
        ec = ec_typer.ECTyper(
            metadata=self.runmetadata,
            report_path=self.report_path,
            assembly_path=os.path.join(self.path, 'raw_assemblies'),
            threads=self.cpus,
            logfile=self.logfile
        )
        ec.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def serosippr(self):
        """
        Serotyping analyses
        """
        Serotype(
            args=self,
            pipelinecommit=self.commit,
            startingtime=self.start_time,
            scriptpath=self.home_path,
            analysistype='serosippr',
            cutoff=0.90,
            pipeline=True
        )
        metadataprinter.MetadataPrinter(inputobject=self)

    def seqsero(self):
        """
        Run SeqSero2 on Salmonella samples
        """
        seqsero_obj = seqsero.SeqSero(self)
        seqsero_obj.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def identity_vtyper(self):
        """
        Legacy vtyper - uses ePCR
        """
        legacy_vtyper = IdentityVtyper(
            metadataobject=self.runmetadata.samples,
            analysistype='legacy_vtyper',
            reportpath=self.report_path,
            mismatches=3
        )
        legacy_vtyper.vtyper()
        metadataprinter.MetadataPrinter(inputobject=self)

    def verotoxin(self):
        """
        Raw read verotoxin typing
        """
        verotoxin = Verotoxin(
            args=self,
            pipeline=True,
            analysistype='verotoxin',
            cutoff=90
        )
        verotoxin.main()

    def sistr(self):
        """
        Sistr
        """
        sistr_obj = sistr.Sistr(
            inputobject=self,
            analysistype='sistr'
        )
        sistr_obj.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def run_gdcs(self):
        """
        Determine the presence of genomically-dispersed conserved sequences
        (genes from MLST, rMLST, and cgMLST analyses)
        """
        # Run the GDCS analysis
        gdcs = GDCS(inputobject=self)
        gdcs.main()
        metadataprinter.MetadataPrinter(inputobject=self)

    def run_report(self):
        """
        Create the final combinedMetadata report
        """
        run_report = reporter.Reporter(self)
        # Create the standard and legacy reports
        run_report.metadata_reporter()
        run_report.legacy_reporter()
        run_report.lab_report()
        # Clean the large attributes from the metadata objects
        run_report.clean_object()

    def __init__(self, args):
        """
        Initialises the variables required for this class
        :param args: list of arguments passed to the script
        """
        self.debug = args.debug
        SetupLogging(self.debug)
        logging.info(
            'Welcome to the CFIA OLC Workflow for Bacterial Assembly and '
            'Typing (COWBAT) version %s', __version__
        )
        # Define variables from the arguments - there may be a more
        # streamlined way to do this
        self.args = args
        if args.sequence_path.startswith('~'):
            self.path = os.path.abspath(
                os.path.expanduser(
                    os.path.join(args.sequence_path)
                )
            )
        else:
            self.path = os.path.abspath(os.path.join(args.sequence_path))
        self.sequence_path = self.path
        if args.reference_file_path.startswith('~'):
            self.ref_file_path = os.path.expanduser(
                os.path.abspath(
                    os.path.join(args.reference_file_path)
                )
            )
        else:
            self.ref_file_path = os.path.abspath(
                os.path.join(args.reference_file_path)
            )
        self.num_reads = args.num_reads
        self.preprocess = args.preprocess
        # Define the start time
        self.start_time = args.startingtime
        if args.custom_sample_sheet:
            if args.custom_sample_sheet.startswith('~'):
                self.custom_sample_sheet = os.path.expanduser(
                    os.path.abspath(
                        os.path.join(self.custom_sample_sheet)
                    )
                )
            else:
                self.custom_sample_sheet = os.path.abspath(
                    os.path.join(args.custom_sample_sheet)
                )
        else:
            self.custom_sample_sheet = args.custom_sample_sheet
        if self.custom_sample_sheet:
            assert os.path.isfile(self.custom_sample_sheet), \
                'Cannot find custom sample sheet as specified ' \
                '{self.custom_sample_sheet}'
        self.basic_assembly = args.basic_assembly
        if not self.custom_sample_sheet and not os.path.isfile(
                os.path.join(self.path, 'SampleSheet.csv')):
            self.basic_assembly = True
            logging.warning(
                'Could not find a sample sheet. Performing basic assembly '
                '(no run metadata captured)'
            )
        # Use the argument for the number of threads to use, or default to the
        # number of cpus in the system
        self.cpus = args.threads if args.threads else \
            multiprocessing.cpu_count() - 1
        # Assertions to ensure that the provided variables are valid
        make_path(self.path)
        assert os.path.isdir(self.path), \
            f'Supplied path location is not a valid directory {self.path!r}'
        self.report_path = os.path.join(self.path, 'reports')
        make_path(self.report_path)
        assert os.path.isdir(self.ref_file_path), \
            f'Reference file path is not a valid directory ' \
            f'{self.ref_file_path!r}'
        self.commit = __version__
        self.home_path = args.home_path
        self.logfile = os.path.join(self.path, 'logfile')
        self.run_info = str()
        self.pipeline = True
        self.qualityobject = MetadataObject()
        # Initialise the metadata object
        self.runmetadata = MetadataObject()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Extract the path of the current script from the full path + file name
    home_path = os.path.split(os.path.abspath(__file__))[0]
    # Parser for arguments
    parser = ArgumentParser(
        description='Assemble genomes from Illumina fastq files'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s commit {__version__}'
    )
    parser.add_argument(
        '-s', '--sequence_path',
        required=True,
        help='Path to folder containing sequencing reads'
    )
    parser.add_argument(
        '-r', '--reference_file_path',
        required=True,
        help='Provide the location of the folder containing the pipeline '
        'accessory files (reference genomes, MLST data, etc.'
    )
    parser.add_argument(
        '-n', '--num_reads',
        default=2,
        type=int,
        help='Specify the number of reads. Paired-reads: 2, unpaired-reads: 1 '
        ' Default is paired-end'
    )
    parser.add_argument(
        '-t', '--threads',
        help='Number of threads. Default is the number of cores in the system'
    )
    parser.add_argument(
        '-c', '--custom_sample_sheet',
        help='Path of folder containing a custom sample sheet and name of '
        'sample sheet file  e.g. /home/name/folder/BackupSampleSheet.csv. '
        'Note that this sheet must still have the same format of Illumina '
        'SampleSheet.csv files'
    )
    parser.add_argument(
        '-b', '--basic_assembly',
        action='store_true',
        help='Performs a basic de novo assembly, and does not collect run '
        'metadata'
    )
    parser.add_argument(
        '-p', '--preprocess',
        action='store_true',
        help='Performs quality trimming and error correction only. Does not '
        'assemble the trimmed + corrected reads'
    )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debug mode for the pipeline (skip FASTQ validation, and '
        'file deletion. Enable debug-level messages if they exist)'
    )
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.startingtime = time()
    arguments.home_path = home_path
    # Run the pipeline
    pipeline = RunAssemble(arguments)
    pipeline.main()
    logging.info('Assembly and characterisation complete')
