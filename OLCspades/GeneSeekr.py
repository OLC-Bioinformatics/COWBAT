#!/usr/bin/env python
import time
import shlex
import subprocess
from csv import DictReader
import copy_reg
import types
from glob import glob
from collections import defaultdict
from Bio.Blast.Applications import NcbiblastnCommandline
from threading import Thread
from accessoryFunctions import *

__author__ = 'mike knowles, adamkoziol'

__doc__ = 'The purpose of this set of modules is to improve upon earlier development of ARMISeekr.py and eventually' \
          'to include generalized functionality for with OOP for GeneSeekr'


# def _pickle_method(method):
#     func_name = method.im_func.__name__
#     obj = method.im_self
#     cls = method.im_class
#     if func_name.startswith('__') and not func_name.endswith('__'):  # deal with mangled names
#         cls_name = cls.__name__.lstrip('_')
#         func_name = '_' + cls_name + func_name
#     return _unpickle_method, (func_name, obj, cls)
#
#
# def _unpickle_method(func_name, obj, cls):
#     for cls in cls.__mro__:
#         try:
#             func = cls.__dict__[func_name]
#         except KeyError:
#             pass
#         else:
#             break
#     return func.__get__(obj, cls)
#
# copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

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

    def makedbthreads(self):
        """
        Setup and create threads for class
        """
        # Find all the target folders in the analysis and add them to the targetfolders set
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
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
            if sample.general.bestassemblyfile != 'NA':
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
                blastn()
            # Run the blast parsing module
            self.blastparser(report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            percentidentity = (float(row['positives']) - float(row['gaps'])) / float(row['subject_length']) * 100
            # Find the allele number and the text before the number for different formats
            target = row['subject_id']
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                self.plusdict[sample.name][target] = percentidentity

    # def makeblastdb(self, (fasta, db)):
    #     print fasta, db
    #     if not os.path.isfile('{}.nhr'.format(db)):  # add nhr for searching
    #         assert os.path.isfile(fasta)  # check that the fasta has been specified properly
    #         MakeBlastDB(db=fasta, out=db, dbtype='nucl')()  # Use MakeBlastDB above
    #     return 0

    # def __init__(self, subject, query, threads=12):
    #     """:type subject: list of genes
    #        :type query: list of target genomes"""
    #     assert isinstance(subject, list), 'Subject is not a list "{0!r:s}"'.format(str(subject))
    #     assert isinstance(query, list), 'Query is not a list"{0!r:s}"'.format(query)
    #     self.count, self.subject, self.query, self.threads = 0, subject, query, threads
    #     self.cutoff, self.genelist = 70, []
    #     self.db = map((lambda x: os.path.splitext(x)[0]), subject)  # remove the file extension for easier globing
    #     self.plus = dict((target, defaultdict(list)) for target in self.query)  # Initialize :return dict
    #     print '[{}] GeneSeekr input is path with {} files'.format(time.strftime("%H:%M:%S"), len(query))
    #     print "[{}] Creating necessary databases for BLAST".format(time.strftime("%H:%M:%S"))
    #     Pool(self.threads).map(self.makeblastdb, zip(self.subject, self.db))
    #     print "\r[{0}] BLAST database(s) created".format(time.strftime("%H:%M:%S"))

    # def _blast(self, (fasta, db)):
    #     blastn = NcbiblastnCommandline(query=fasta,
    #                                    db=db,
    #                                    evalue=10,
    #                                    outfmt="'6 sseqid nident slen'",
    #                                    perc_identity=self.cutoff)
    #     stdout, stderr = blastn()
    #     if stdout != '':
    #         return [[fasta, aln[0], '{:.2f}'.format(abs(float(aln[1]) / float(aln[2]) * 100))]
    #                 for aln in [hsp.split('\t')
    #                             for hsp in stdout.rstrip().split("\n")]
    #                 if abs(float(aln[1]) / float(aln[2])) >= self.cutoff/100.0]
    #
    # def mpblast(self, cutoff=70):
    #     assert isinstance(cutoff, int), 'Cutoff is not an integer {0!r:s}'.format(str(cutoff))
    #     self.cutoff = cutoff
    #     print "[{}] Now performing and parsing BLAST database searches".format(time.strftime("%H:%M:%S"))
    #     start = time.time()
    #     p = Pool(12)
    #     for genes in self.db:
    #         mapblast = p.map(self._blast, [(genome, genes) for genome in self.query])
    #         for fastaline in mapblast:
    #             if fastaline is not None:  # if the returned list contains [genome, gene, value]
    #                 for fasta, gene, v in fastaline:  # unpack
    #                     if gene not in self.genelist:
    #                         self.genelist.append(gene)  # create list of all genes in anaylsis
    #                     self.plus[fasta][gene].append(v)
    #                     self.plus[fasta][gene].sort()
    #
    #     end = time.time() - start
    #     print "\n[{0:s}] Elapsed time for GeneSeekr is {1:0d}m {2:0d}s with {3:0.2f}s per genome".format(
    #         time.strftime("%H:%M:%S"), int(end) / 60, int(end) % 60, end / float(len(self.query)))
    #
    # def csvwriter(self, out, name):
    #     assert isinstance(out, str), u'Output location is not a string {0!r:s}'.format(out)
    #     assert isinstance(name, str), u'Output name is not a string {0!r:s}'.format(name)
    #     assert os.path.isdir(out), 'Output location is not a valid directory {0!r:s}'.format(out)
    #     self.genelist.sort()
    #     rowcount, row = 0, 'Strain,'
    #     row += ', '.join(self.genelist)
    #     for genomerow in sorted(self.plus):
    #         row += '\n{}'.format(os.path.split(os.path.splitext(genomerow)[0])[1].replace('_filteredAssembled', ""))
    #         for genename in self.genelist:
    #             row += ',' + (lambda x, y: ' '.join(map(str, x[y])) if y in x else 'N')(self.plus[genomerow],
    # genename)
    #             # Add the allele numbers to the row for the appropriate gene, otherwise return N
    #     with open("%s/%s_results_%s.csv" % (out, name, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
    #         csvfile.write(row)

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = inputobject.start
        self.analysistype = inputobject.analysistype
        self.targetfolders = set()
        # Fields used for custom outfmt 6 BLAST output:
        # "6 qseqid sseqid positive mismatch gaps evalue bitscore slen length"
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue',  'bit_score', 'subject_length', 'alignment_length']
        #
        self.plusdict = defaultdict(make_dict)
        self.dqueue = Queue()
        self.blastqueue = Queue()
        self.geneseekr()


# def helper(genes, targets, out, cuttoff, threads):
#     from glob import glob
#     assert os.path.isdir(out), u'Output location is not a valid directory {0!r:s}'.format(out)
#     assert os.path.isfile(genes), u'ARMI-genes.fa not valid {0!r:s}'.format(genes)
#     # assert os.path.isfile(aro), u'Antibiotic JSON not valid {0!r:s}'.format(aro)
#     assert isinstance(threads, int)
#     ispath = (lambda x: glob(x + "/*.fasta") if os.path.isdir(x) else [x])
#     genes = ispath(genes)
#     targets = ispath(targets)
#     result = GeneSeekr(genes, targets, threads)
#     result.mpblast(cuttoff)
#     result.csvwriter(out, 'test')
    # json.dump(result.plus, open("%s/ARMI-gene_results_%s.json" % (out, time.strftime("%Y.%m.%d.%H.%M.%S")), 'w'),
    #           sort_keys=True, indent=4, separators=(',', ': '))
    # decipher(result.plus, json.load(open(aro)), out)


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
            parser = ArgumentParser(description='Antibiotic Resistance Marker Identifier:\n'
                                                'Use to find markers for any bacterial genome')
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
            # args['anti'],
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
            # Run the analyses
            GeneSeekr(self)

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                targetpath = '{}{}/{}'.format(self.referencefilepath, self.analysistype, sample.general.referencegenus)
                targets = glob('{}/*.fa*'.format(targetpath))
                try:
                    combinedtargets = glob('{}/*.tfa'.format(targetpath))[0]
                except IndexError:
                    combinetargets(targets, targetpath)
                    combinedtargets = glob('{}/*.tfa'.format(targetpath))
                sample[self.analysistype].targets = targets
                sample[self.analysistype].combinedtargets = combinedtargets
                sample[self.analysistype].targetpath = targetpath
                sample[self.analysistype].targetnames = [os.path.split(x)[1].split('.')[0] for x in targets]
                sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory, self.analysistype)
            else:
                # Set the metadata file appropriately
                sample[self.analysistype].targets = 'NA'
                sample[self.analysistype].combinedtargets = 'NA'
                sample[self.analysistype].targetpath = 'NA'
                sample[self.analysistype].targetnames = 'NA'
                sample[self.analysistype].reportdir = 'NA'

    def __init__(self, inputobject, analysistype):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.referencefilepath = inputobject.reffilepath
        self.reportdir = '{}/'.format(inputobject.reportpath)
        self.cutoff = 70
        # Get the alleles and profile into the metadata
        self.strainer()
        GeneSeekr(self)

# TODO no pluses
# TODO multifasta?
# TODO alignments
