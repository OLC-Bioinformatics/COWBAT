#!/usr/bin/env python
import copyreg
import types
import os
from GeneSeekr import *
from accessoryFunctions import make_path

__author__ = 'adamkoziol'


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):  # deal with mangled names
        cls_name = cls.__name__.lstrip('_')
        func_name = '_' + cls_name + func_name
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.__mro__:
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    # noinspection PyUnboundLocalVariable
    return func.__get__(obj, cls)

copyreg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class ARMI(GeneSeekr):

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            size = 0
            make_path(sample[self.analysistype].reportdir)
            try:
                report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
                size = os.path.getsize(report)
            except IndexError:

                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
            db = target.split('.')[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly, db=db,
                                           num_alignments=1000000,
                                           num_threads=12,
                                           evalue=1e-4,
                                           outfmt="'6 sseqid nident gaps slen qacc qstart qend'",
                                           out=report)
            # Save the blast command in the metadata
            sample[self.analysistype].blastcommand = str(blastn)
            if not os.path.isfile(report) or size == 0:
                blastn()
            # Run the blast parsing module
            self.blastparser(report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.armifields, dialect='excel-tab')
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = abs(float('{:0.2f}'.format((float(row['num_identities']) - float(row['gaps'])) /
                                                         float(row['subject_length']) * 100)))
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Pull the gene and allele from the subject_id
                gene = row['subject_id'].split('|')[0][4:]
                allele = row['subject_id'].split('|')[1]
                # Don't add duplicate hits to the dictionary
                if gene not in self.plus[sample.name]:
                    # Append the sample name, the first seven characters of the gene, and a list of the percent
                    # identity, the query accession, and the allele name to the dictionary
                    self.plus[sample.name][gene[:7]].append([percentidentity, row['query_accession'], allele])

    def csvwriter(self):
        from .ARMICARD import decipher
        import pickle
        aro = ''
        # Set the location of the .dat file - it is inefficiently done for each sample, as I don't have the .targetpath
        # stored anywhere but the metadata
        for sample in self.metadata:
            if sample[self.analysistype].targetpath != 'NA' and sample.general.bestassemblyfile != 'NA'\
                    and sample[self.analysistype].reportdir != 'NA':
                aro = '{}aro.dat'.format(sample[self.analysistype].targetpath)
        # Send the dictionaries and variables to decipher - it will analyse the BLAST results and create reports
        if os.path.isfile(aro):
            with open(aro) as anti:
                decipher(self.plus, pickle.load(anti), self.reportpath, self.metadata, tolc=False)

    def __init__(self, inputobject):
        # '6 sseqid nident gaps slen qacc qstart qend'

        self.armifields = ['subject_id', 'num_identities', 'gaps', 'subject_length', 'query_accession', 'query_start',
                           'query_end']
        # Initialise dictionary with the sample names
        self.plus = dict((sample.name, defaultdict(list)) for sample in inputobject.runmetadata.samples)
        # Run the GeneSeekr __init__
        GeneSeekr.__init__(self, inputobject)
