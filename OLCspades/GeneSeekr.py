#!/usr/bin/env python
import time
import shlex
import subprocess
from csv import DictReader
from glob import glob
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline
from threading import Thread
from accessoryFunctions import *

__author__ = 'mike knowles, adamkoziol'

__doc__ = 'The purpose of this set of modules is to improve upon earlier development of ARMISeekr.py and eventually' \
          'to include generalized functionality for with OOP for GeneSeekr'


# from mMLST import MLST
class GeneSeekr(object):

    def geneseekr(self):
        # Make blast databases (if necessary)
        printtime('Creating {} blast databases as required'.format(self.analysistype), self.start)
        self.makedbthreads()
        # Run the blast analyses
        printtime('Running {} blast analyses'.format(self.analysistype), self.start)
        self.blastnthreads()
        globalcounter()
        self.csvwriter()
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata:
            delattr(sample[self.analysistype], "targetnames")
            delattr(sample[self.analysistype], "targets")
        printtime('{} analyses complete'.format(self.analysistype), self.start)

    def makedbthreads(self):
        """
        Setup and create threads for class
        """
        # Find all the target folders in the analysis and add them to the targetfolders set
        for sample in self.metadata:
            if sample[self.analysistype].combinedtargets != 'NA':
                self.targetfolders.add(sample[self.analysistype].targetpath)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.targetfolders)):
            # Send the threads to makeblastdb
            threads = Thread(target=self.makeblastdb, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Make blast databases for MLST files (if necessary)
        for targetdir in self.targetfolders:
            # List comprehension to remove any previously created database files from list
            targetfiles = glob('{}/*.tfa'.format(targetdir))
            for targetfile in targetfiles:
                # Add the fasta file to the queue
                self.dqueue.put(targetfile)
        self.dqueue.join()  # wait on the dqueue until everything has been processed

    def makeblastdb(self):
        """Makes blast database files from targets as necessary"""
        while True:  # while daemon
            fastapath = self.dqueue.get()  # grabs fastapath from dqueue
            # remove the path and the file extension for easier future globbing
            db = fastapath.split('.')[0]
            nhr = '{}.nhr'.format(db)  # add nhr for searching
            fnull = open(os.devnull, 'w')  # define /dev/null
            if not os.path.isfile(str(nhr)):  # if check for already existing dbs
                # Create the databases
                # TODO use MakeBLASTdb class
                subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                            .format(fastapath, db)), stdout=fnull, stderr=fnull)
                # os.system('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                #           .format(fastapath, db))
            dotter()
            self.dqueue.task_done()  # signals to dqueue job is done

    def blastnthreads(self):
        """Setup and create  threads for blastn and xml path"""
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                for i in range(len(sample[self.analysistype].combinedtargets)):
                    threads = Thread(target=self.runblast, args=())
                    threads.setDaemon(True)
                    threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample[self.analysistype].combinedtargets != 'NA':
                # Add each fasta file combination to the threads
                self.blastqueue.put((sample.general.bestassemblyfile, sample[self.analysistype].combinedtargets,
                                     sample))
        # Join the threads
        self.blastqueue.join()

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            make_path(sample[self.analysistype].reportdir)
            try:
                report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
                size = os.path.getsize(report)
                if size == 0:
                    os.remove(report)
                    report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                              time.strftime("%Y.%m.%d.%H.%M.%S"))
            except IndexError:
                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
            db = target.split('.')[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly, db=db, evalue='1E-20', num_alignments=1000000,
                                           num_threads=12,
                                           outfmt='"6 qseqid sseqid positive mismatch gaps '
                                                  'evalue bitscore slen length"',
                                           out=report)
            if not os.path.isfile(report):
                try:
                    blastn()
                except:
                    self.blastqueue.task_done()
                    self.blastqueue.join()
                    try:
                        os.remove(report)
                    except IOError:
                        pass
                    raise
            # Run the blast parsing module
            self.blastparser(report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        resultdict = {}
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            target = row['subject_id']
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Update the dictionary with the target and percent identity
                resultdict.update({target: percentidentity})
            sample[self.analysistype].blastresults = resultdict
        if not resultdict:
            sample[self.analysistype].blastresults = 'NA'

    def csvwriter(self):
        combinedrow = ''
        for sample in self.metadata:
            row = ''
            if sample[self.analysistype].targetnames != 'NA':
                # Populate the header with the appropriate data, including all the genes in the list of targets
                row += 'Strain,{},\n'.format(','.join(sorted(sample[self.analysistype].targetnames)))
                row += '{},'.format(sample.name)
                for target in sorted(sample[self.analysistype].targetnames):
                    if sample[self.analysistype].blastresults != 'NA':
                        try:
                            type(sample[self.analysistype].blastresults[target])
                            row += '{},'.format(sample[self.analysistype].blastresults[target])
                        except (KeyError, TypeError):
                            row += '-,'
                    else:
                        row += '-,'
                row += '\n'
                combinedrow += row
                # If the length of the number of report directories is greater than 1 (script is being run as part of
                # the assembly pipeline) make a report for each sample
                if self.pipeline:
                    # Open the report
                    with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                   self.analysistype), 'wb') as report:
                        # Write the row to the report
                        report.write(row)
            else:
                sample[self.analysistype].blastresults = 'NA'

        # Create the report containing all the data from all samples
        if self.pipeline:
            with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                    as combinedreport:
                combinedreport.write(combinedrow)
        else:
            with open('{}{}_{:}.csv'.format(self.reportpath, self.analysistype, time.strftime("%Y.%m.%d.%H.%M.%S")),
                      'wb') \
                    as combinedreport:
                combinedreport.write(combinedrow)

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = inputobject.start
        self.analysistype = inputobject.analysistype
        self.reportpath = inputobject.reportdir
        self.targetfolders = set()
        self.pipeline = inputobject.pipeline
        self.referencefilepath = inputobject.referencefilepath
        # Fields used for custom outfmt 6 BLAST output:
        # "6 qseqid sseqid positive mismatch gaps evalue bitscore slen length"
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue',  'bit_score', 'subject_length', 'alignment_length']
        #
        self.plusdict = defaultdict(make_dict)
        self.dqueue = Queue()
        self.blastqueue = Queue()
        self.geneseekr()


if __name__ == '__main__':

    class Parser(object):

        def strainer(self):
            from accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = sorted(glob('{}*.fa*'.format(self.sequencepath)))
            self.targets = sorted(glob('{}*.fa*'.format(self.targetpath)))
            try:
                self.combinedtargets = glob('{}/*.tfa'.format(self.targetpath))[0]
            except IndexError:
                combinetargets(self.targets, self.targetpath)
                self.combinedtargets = glob('{}/*.tfa'.format(self.targetpath))[0]
            # Populate the metadata object. This object will be populated to mirror the objects created in the
            # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
            # of a pipeline
            assert self.strains, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your sequence path is correct'.format(self.sequencepath)
            assert self.targets, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your target path is correct'.format(self.targetpath)
            for sample in self.strains:
                # Create the object
                metadata = MetadataObject()
                # Set the base file name of the sequence. Just remove the file extension
                filename = os.path.split(sample)[1].split('.')[0]
                # Set the .name attribute to be the file name
                metadata.name = filename
                # Create the .general attribute
                metadata.general = GenObject()
                # Create the .mlst attribute
                setattr(metadata, self.analysistype, GenObject())
                # metadata.mlst = GenObject()
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = sample
                metadata[self.analysistype].targets = self.targets
                metadata[self.analysistype].combinedtargets = self.combinedtargets
                metadata[self.analysistype].targetpath = self.targetpath
                metadata[self.analysistype].targetnames = [os.path.split(x)[1].split('.')[0] for x in self.targets]
                metadata[self.analysistype].reportdir = self.reportpath
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            parser = ArgumentParser(description='Use to find markers for any bacterial genome')
            parser.add_argument('--version',
                                action='version',
                                version='%(prog)s v0.5')
            parser.add_argument('-s', '--sequencepath',
                                required=True,
                                help='Specify input fasta folder')
            parser.add_argument('-t', '--targetpath',
                                required=True,
                                help='Specify folder of targets')
            parser.add_argument('-r', '--reportpath',
                                required=True,
                                help='Specify output folder for csv')
            # parser.add_argument('-a', '--anti', type=str, required=True, help='JSON file location')
            parser.add_argument('-c', '--cutoff',
                                type=int,
                                default=70, help='Threshold for maximum unique bacteria for a single antibiotic')
            parser.add_argument('-n', '--numthreads',
                                type=int,
                                default=24,
                                help='Specify number of threads')
            args = parser.parse_args()
            self.sequencepath = os.path.join(args.sequencepath, '')
            assert os.path.isdir(self.sequencepath), 'Cannot locate sequence path as specified: {}'\
                .format(self.sequencepath)
            self.targetpath = os.path.join(args.targetpath, '')
            assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}'\
                .format(self.targetpath)
            self.reportpath = os.path.join(args.reportpath, '')
            assert os.path.isdir(self.reportpath), 'Cannot locate report path as specified: {}'\
                .format(self.reportpath)
            self.cutoff = args.cutoff
            self.threads = args.numthreads
            #
            self.strains = []
            self.targets = []
            self.combinedtargets = ''
            self.samples = []
            self.analysistype = 'geneseekr'
            self.start = time.time()
            self.strainer()

    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.start = self.runmetadata.start
            self.analysistype = self.runmetadata.analysistype
            # self.alleles = self.runmetadata.alleles
            # self.profile = self.runmetadata.profile
            self.cutoff = self.runmetadata.cutoff
            self.threads = self.runmetadata.threads
            self.reportdir = self.runmetadata.reportpath
            self.pipeline = False
            self.referencefilepath = ''
            # Run the analyses
            GeneSeekr(self)

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                if self.genusspecific:
                    targetpath = '{}{}/{}'.format(self.referencefilepath, self.analysistype,
                                                  sample.general.referencegenus)
                else:
                    targetpath = '{}{}/'.format(self.referencefilepath, self.analysistype)
                targets = glob('{}/*.fa*'.format(targetpath))
                targetcheck = glob('{}/*.*fa*'.format(targetpath))
                if targetcheck:
                    try:
                        combinedtargets = glob('{}/*.tfa'.format(targetpath))[0]
                    except IndexError:
                        combinetargets(targets, targetpath)
                        combinedtargets = glob('{}/*.tfa'.format(targetpath))[0]
                    sample[self.analysistype].targets = targets
                    sample[self.analysistype].combinedtargets = combinedtargets
                    sample[self.analysistype].targetpath = targetpath
                    sample[self.analysistype].targetnames = [os.path.split(x)[1].split('.')[0] for x in targets]
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                else:
                    # Set the metadata file appropriately
                    sample[self.analysistype].targets = 'NA'
                    sample[self.analysistype].combinedtargets = 'NA'
                    sample[self.analysistype].targetpath = 'NA'
                    sample[self.analysistype].targetnames = 'NA'
                    sample[self.analysistype].reportdir = 'NA'
            else:
                # Set the metadata file appropriately
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].targets = 'NA'
                sample[self.analysistype].combinedtargets = 'NA'
                sample[self.analysistype].targetpath = 'NA'
                sample[self.analysistype].targetnames = 'NA'
                sample[self.analysistype].reportdir = 'NA'

    def __init__(self, inputobject, analysistype, genusspecific, cutoff):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.referencefilepath = inputobject.reffilepath
        self.reportdir = '{}/'.format(inputobject.reportpath)
        self.cutoff = cutoff
        self.pipeline = True
        self.genusspecific = genusspecific
        # Get the alleles and profile into the metadata
        self.strainer()
        # GeneSeekr(self)

# TODO no pluses
# TODO multifasta?
# TODO alignments
