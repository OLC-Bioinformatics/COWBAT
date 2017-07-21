#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class Serotype(GeneSeekr):

    def reporter(self):
        combinedrow = ''
        for sample in self.metadata:
            row = ''
            try:
                # Serotyping schemes differ depending on the organism
                if sample.general.referencegenus == 'Escherichia':
                    row += 'Strain,O-type,H-type\n'
                    # Make lists of all the tuples in the dictionary with queries that, following splitting on
                    # underscore have their final entry begin with O, or with H, respectively (e.g. fliC_56_AJ884569_H8
                    # would yield H8). Sort the lists based on the second entry in each tuple
                    # (reverse makes the list descending)
                    otypes = sorted([item for item in sample[self.analysistype].blastresults.items()
                                     if list(item)[0].split('_')[-1].startswith('O')], key=lambda x: x[1], reverse=True)
                    htypes = sorted([item for item in sample[self.analysistype].blastresults.items()
                                     if list(item)[0].split('_')[-1].startswith('H')], key=lambda x: x[1], reverse=True)
                    # Try/except loops to set the O/H-type and the identity for that O/H-type for each sample
                    try:
                        otype = otypes[0][0].split('_')[-1]
                        oidentity = '({})'.format(otypes[0][1])
                    except IndexError:
                        otype = '-'
                        oidentity = ''
                    try:
                        htype = htypes[0][0].split('_')[-1]
                        hidentity = '({})'.format(htypes[0][1])
                    except IndexError:
                        htype = '-'
                        hidentity = ''
                    sample[self.analysistype].serotype = '{}:{}'.format(otype, htype)
                    # Populate the row for the report
                    row += '{},{}{},{}{}\n'.format(sample.name, otype, oidentity, htype, hidentity)
                # Will be enhanced when more serotyping schemes are incorporated
                else:
                    sample[self.analysistype].serotype = 'NA'
                combinedrow += row
                # If the script is being run as part of the assembly pipeline, make a report for each sample
                if self.pipeline:
                    if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].reportdir != 'NA':
                        # Open the report
                        with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                       self.analysistype), 'w') as report:
                            # Write the row to the report
                            report.write(row)
                # Remove the messy blast results from the object
                delattr(sample[self.analysistype], "blastresults")
            except (AttributeError, KeyError):
                sample[self.analysistype].serotype = 'NA'
                pass
        # Create the report containing all the data from all samples
        with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'w') \
                as combinedreport:
            combinedreport.write(combinedrow)
