#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import MetadataObject, GenObject, make_path, relative_symlink, SetupLogging
from spadespipeline.legacy_vtyper import Vtyper as LegacyVtyper
import accessoryFunctions.metadataprinter as metadataprinter
from spadespipeline.typingclasses import Prophages, Univec
from spadespipeline.createobject import ObjectCreation
from spadespipeline.mobrecon import MobRecon
from metagenomefilter import automateCLARK
import spadespipeline.prodigal as prodigal
import spadespipeline.reporter as reporter
import spadespipeline.quality as quality
import spadespipeline.sistr as sistr
from geneseekr.blast import BLAST
import coreGenome.core as core
import MASHsippr.mash as mash
from argparse import ArgumentParser
import multiprocessing
from time import time
import logging
import os

__version__ = '0.0.01'
__author__ = 'adamkoziol'


class Typing(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        # Create the metadata objects
        self.objects()
        # Determine assembly stats
        self.assembly_stats()
        # Perform genus-agnostic typing
        self.agnostictyping()
        # Perform typing
        self.typing()
        #
        self.typing_reports()
        # Compress or remove all large, temporary files created by the pipeline
        # compress.Compress(self)
        metadataprinter.MetadataPrinter(inputobject=self)

    def objects(self):
        """

        :return:
        """
        self.runmetadata = ObjectCreation(inputobject=self)
        make_path(os.path.join(self.path, 'BestAssemblies'))
        for sample in self.runmetadata.samples:
            # Link the assemblies to the BestAssemblies folder - necessary for GenomeQAML
            relative_symlink(sample.general.bestassemblyfile,
                             os.path.join(self.path, 'BestAssemblies'))
            # Create attributes required for downstream analyses
            sample.general.trimmedcorrectedfastqfiles = [sample.general.bestassemblyfile]
        # self.metadata = self.metadata.samples

    def assembly_stats(self):
        """
        Perform some basic quality analyses on the assemblies
        """
        # Calculate assembly metrics on polished assemblies
        self.quality_features(analysis='polished')
        # ORF detection
        self.prodigal()
        # Assembly quality determination
        self.genome_qaml()
        # CLARK analyses
        self.clark()
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
        self.geneseekr()
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
        logging.info('Running MASH analyses')
        mash.Mash(inputobject=self,
                  analysistype='mash')
        metadataprinter.MetadataPrinter(inputobject=self)

    def rmlst(self):
        """
        Run rMLST analyses
        """
        logging.info('rMLST')
        rmlst = BLAST(args=self,
                      analysistype='rMLST',
                      cutoff=100,
                      unique=True,
                      pipeline=True)
        rmlst.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def sixteens(self):
        """
        Run the 16S analyses
        """
        sixteen_s = BLAST(args=self,
                          analysistype='sixteens_full',
                          cutoff=95)
        sixteen_s.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def run_gdcs(self):
        """
        Determine the presence of genomically-dispersed conserved sequences for Escherichia, Listeria, and Salmonella
        strains
        """
        # Run the GDCS analysis
        gdcs = BLAST(args=self,
                     analysistype='GDCS',
                     cutoff=90,
                     genus_specific=True,
                     unique=False,)
        gdcs.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def geneseekr(self):
        """
        Find genes of interest
        """
        geneseekr = BLAST(args=self,
                          analysistype='genesippr',
                          cutoff=95)
        geneseekr.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def mob_suite(self):
        """

        """
        mob = MobRecon(metadata=self.runmetadata.samples,
                       analysistype='mobrecon',
                       databasepath=self.targetpath,
                       threads=self.cpus,
                       logfile=self.logfile,
                       reportpath=self.reportpath)
        mob.mob_recon()
        metadataprinter.MetadataPrinter(inputobject=self)

    def resfinder(self):
        """
        Resistance finding - assemblies
        """
        resfinder = BLAST(args=self,
                          analysistype='resfinder_assembled')
        resfinder.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def prophages(self):
        """
        Prophage detection
        """
        prophages = Prophages(args=self,
                              analysistype='prophages',
                              cutoff=90,
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
        virulence = BLAST(args=self,
                          analysistype='virulence')
        virulence.seekr()
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
        logging.info('MLST')
        mlst = BLAST(args=self,
                     analysistype='MLST',
                     genus_specific=True,
                     cutoff=100,
                     unique=True)
        mlst.seekr()
        metadataprinter.MetadataPrinter(inputobject=self)

    def serosippr(self):
        """
        Serotyping analyses
        """
        #          pipeline=True)
        sero = BLAST(args=self,
                     analysistype='serosippr',
                     cutoff=90,
                     genus_specific=True,
                     unique=True)
        sero.seekr()
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

    def typing_reports(self):
        """
        Create empty attributes for analyses that were not performed, so that the metadata report can be created
        :return:
        """
        for sample in self.runmetadata.samples:
            sample.confindr = GenObject()
            sample.mapping = GenObject()
            sample.mapping.MeanInsertSize = 'ND'
            sample.mapping.MeanCoveragedata = 'ND'
            sample.mlst.results = sample.mlst.blastresults
            sample.genesippr.report_output = set()
            sample.genesippr.results = dict()
            try:
                for gene, percentid in sample.genesippr.blastresults.items():
                    if percentid > 95:
                        sample.genesippr.report_output.add(gene.split('_')[0])
            except AttributeError:
                sample.genesippr.report_output = list()
            sample.genesippr.report_output = sorted(list(sample.genesippr.report_output))
        # Create a report
        reporter.Reporter(inputobject=self,
                          legacy=False)

    def __init__(self, start, sequencepath, referencefilepath, scriptpath):
        """
        
        :param start: 
        :param sequencepath: 
        :param referencefilepath: 
        :param scriptpath:
        """
        SetupLogging(debug=True)
        logging.info('Welcome to the CFIA bacterial typing pipeline {}'
                     .format(__version__))
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.sequencepath = os.path.join(sequencepath)
        self.path = self.sequencepath
        self.targetpath = os.path.join(referencefilepath)
        self.reffilepath = self.targetpath
        # Define the start time
        self.starttime = start
        self.start = self.starttime
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = multiprocessing.cpu_count() - 1
        # Assertions to ensure that the provided variables are valid
        assert os.path.isdir(self.sequencepath), 'Supplied path location is not a valid directory {0!r:s}'\
            .format(self.sequencepath)
        self.reportpath = os.path.join(self.sequencepath, 'reports')
        assert os.path.isdir(self.targetpath), 'Reference file path is not a valid directory {0!r:s}'\
            .format(self.targetpath)
        self.commit = __version__
        self.homepath = scriptpath
        self.analysistype = 'assembly_typing'
        self.genus_specific = False
        self.logfile = os.path.join(self.sequencepath, 'logfile')
        self.pipeline = True
        # Initialise the metadata object
        self.metadata = list()
        self.runmetadata = MetadataObject()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Parser for arguments
    parser = ArgumentParser(description='Performing the typing component of the COWBAT pipeline on assemblies')
    parser.add_argument('-v', '--version',
                        action='version', version=__version__)
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to folder containing sequencing reads')
    parser.add_argument('-r', '--referencefilepath',
                        required=True,
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    arguments = parser.parse_args()
    # Run the pipeline
    pipeline = Typing(start=time(),
                      sequencepath=arguments.sequencepath,
                      referencefilepath=arguments.referencefilepath,
                      scriptpath=homepath)
    pipeline.main()
    logging.info('Characterisation complete')
