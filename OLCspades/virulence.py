#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class Virulence(GeneSeekr):

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
            # Find the allele number and the text before the number for different formats
            target = '{},{}'.format(row['subject_id'], row['query_id'])
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # self.plusdict[sample.name][target] = percentidentity
                resultdict.update({target: percentidentity})
            sample[self.analysistype].blastresults = resultdict

    def csvwriter(self):
        combinedrow = ''
        for sample in self.metadata:
            row = ''
            # Virulence schemes differ depending on the organism
            if sample.general.referencegenus == 'Escherichia':
                row += 'Contig,Gene,Percentidentity\n'
                # print sample[self.analysistype].blastresults.items()
                geneset = sorted(set([gene[0].split(':')[0] for gene in
                                      sample[self.analysistype].blastresults.items()]))
                # print sample.name, geneset
                for gene in geneset:
                    # print sample.name, gene
                    besthits = sorted([item for item in sample[self.analysistype].blastresults.items()
                                       if list(item)[0].startswith(gene)], key=lambda x: x[1], reverse=True)
                    # print gene, besthits
                    if 'stx' in gene:
                        # print besthits
                        # print besthits[0][0].split(',')[0].split(':')[-1]
                        for hit in besthits:
                            if hit[1] > 98 and '{}{}'.format(gene, hit[0].split(',')[0].split(':')[-1]) not in row:
                                # print hit, hit[0].split(',')[0].split(':')[-1]
                                if len(hit[0].split(',')[0].split(':')[-1]) == 1:
                                    row += '{},{}{},{}\n'\
                                        .format(hit[0].split(',')[1], gene, hit[0].split(',')[0].split(':')[-1], hit[1])
                                elif gene not in row:
                                    row += '{},{},{}\n'.format(hit[0].split(',')[1], gene, hit[1])
                        sample.general.stx = True
                        # print besthits[0][1]
                    else:
                        row += '{},{},{}\n'.format(besthits[0][0].split(',')[1], gene, besthits[0][1])
                        if not sample.general.stx:
                            sample.general.stx = False

            combinedrow += row
            # If the script is being run as part of the assembly pipeline, make a report for each sample
            if self.pipeline:
                # Open the report
                with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                               self.analysistype), 'wb') as report:
                    # Write the row to the report
                    report.write(row)
            # Remove the messy blast results from the object
            delattr(sample[self.analysistype], "blastresults")
        # Create the report containing all the data from all samples
        with open('{}{}_{:}.csv'.format(self.reportpath, self.analysistype, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') \
                as combinedreport:
            combinedreport.write(combinedrow)
