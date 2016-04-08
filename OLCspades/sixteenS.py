#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class SixteenS(GeneSeekr):

    def csvwriter(self):
        import operator
        from glob import glob
        from csv import DictReader
        # The field names of the columns in the taxonomy file
        fieldnames = ['Accession', 'Taxonomy']
        # Open the sequence profile file as a dictionary
        for sample in self.metadata:
            # Sort the results stored in the object and pull out the best match
            try:
                bestmatch = sorted(sample[self.analysistype].blastresults.items(),
                                   key=operator.itemgetter(1), reverse=True)[0]
                # The taxonomy file is a tab-delimited file with two columns: Accession and Taxonomy
                taxafile = glob('{}*.tax'.format(sample[self.analysistype].targetpath))[0]
                # Open the taxonomy file as a dictionary
                taxondict = DictReader(open(taxafile), fieldnames=fieldnames, dialect='excel-tab')
                # The accession supplied in the blast output merges the .1 at the end AF403733.1 is displayed as
                # AF4037331. Take all the characters in the accession except the last, and add a .1 to recreate
                # the proper accession
                fixedaccession = bestmatch[0][:-1] + ".1"
                for row in taxondict:
                    if row['Accession'] == fixedaccession:
                        sample[self.analysistype].taxonomy = row['Taxonomy']
                        sample[self.analysistype].accession = fixedaccession
                        sample[self.analysistype].percentidentity = bestmatch[1]
                sample[self.analysistype].taxafile = taxafile
            except (AttributeError, KeyError):
                sample[self.analysistype].taxonomy = 'NA'
                sample[self.analysistype].accession = 'NA'
                sample[self.analysistype].percentidentity = 'NA'
                sample[self.analysistype].taxafile = 'NA'
            try:
                # Remove the messy blast results from the object
                delattr(sample[self.analysistype], "blastresults")
            except KeyError:
                pass
