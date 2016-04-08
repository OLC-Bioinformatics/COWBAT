#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class CoreGenome(GeneSeekr):

    def csvwriter(self):
        combinedrow = ''
        for sample in self.metadata:
            header = 'Strain,GenesPresent/Total,'
            row = '{},'.format(sample.name)
            data = ''
            # print sample.name, sample[self.analysistype].targetnames
            try:
                count = 0
                for gene in sample[self.analysistype].targetnames:
                    header += '{},'.format(gene)
                # for result in sample[self.analysistype].blastresults.items():
                    try:
                        besthit = sorted([[item.split('|')[0], sample[self.analysistype].blastresults[item]] for item
                                          in sample[self.analysistype].blastresults if gene in item],
                                         key=lambda x: x[1], reverse=True)[0]
                        data += '{},'.format(besthit[1])
                        # print sample.name, gene, besthit
                        count += 1
                        # row +=
                    except IndexError:
                        data += '-,'
                # print count, len(sample[self.analysistype].targetnames)
                header += '\n'
                row += '{}/{},'.format(count, len(sample[self.analysistype].targetnames))
                row += data
                sample[self.analysistype].targetspresent = count
                sample[self.analysistype].totaltargets = len(sample[self.analysistype].targetnames)
                percentage = float(count) / float(len(sample[self.analysistype].targetnames)) * 100
                sample[self.analysistype].targetpercentage = '{:0.2f}'.format(percentage)

                # If the script is being run as part of the assembly pipeline, make a report for each sample
                if self.pipeline:
                    # Open the report
                    with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                   self.analysistype), 'wb') as report:
                        # Write the row to the report
                        report.write(header)
                        report.write(row)
                # Remove the messy blast results from the object
                try:
                    delattr(sample[self.analysistype], "blastresults")
                except KeyError:
                    pass
                    # print sample.name, gene
            except KeyError:
                sample[self.analysistype].targetspresent = 'NA'
                sample[self.analysistype].totaltargets = 'NA'
                sample[self.analysistype].targetpercentage = 'NA'

