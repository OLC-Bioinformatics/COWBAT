#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
__author__ = 'adamkoziol'


class Submit(object):

    def submit(self):
        for sample in self.runmetadata.samples:
            if sample[self.analysistype].reportdir != 'NA':
                #
                if sample[self.analysistype].closealleles:
                    #
                    sample[self.analysistype].closesequences = dict()
                    for gene, allele in sample[self.analysistype].closealleles.items():
                        dnaseq = Seq(sample[self.analysistype].queryseq[gene], IUPAC.unambiguous_dna)
                        protseq = str(dnaseq.translate())
                        print gene, allele, sample[self.analysistype].start[gene], sample[self.analysistype].end[gene],\
                            sample[self.analysistype].mismatches[gene], sample[self.analysistype].querylength[gene], \
                            sample[self.analysistype].subjectlength[gene], sample[self.analysistype].queryid[gene]
                        print sample[self.analysistype].queryseq[gene]
                        print protseq
                        print protseq.count('*')

    def __init__(self, inputobject, analysistype):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.start
        self.referencefilepath = inputobject.referencefilepath
        self.reportdir = inputobject.reportdir
        #
        self.submit()
