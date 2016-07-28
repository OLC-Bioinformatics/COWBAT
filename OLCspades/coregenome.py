#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class CoreGenome(GeneSeekr):

    def blastparser(self, report, sample):
        import operator
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        resultdict = {}
        # Create a list of all the names of the database files - glob, replace - with _, remove path and extension
        coregenomes = map(lambda x: os.path.basename(x).split('.')[0].replace('-', '_'),
                          glob('{}/*.fas'.format(os.path.join(self.referencefilepath,
                                                              self.analysistype, sample.general.referencegenus))))
        # Initialise resultdict with an integer for every database file
        for genome in coregenomes:
            resultdict[genome] = int()
        # A set to store the number of core genes
        coregenes = set()
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            # Split off any | from the sample name
            target = row['subject_id'].split('|')[0]
            # As there are variable numbers of _ in the name, try to split off the last one only if there are multiple
            # and only keep the first part of the split if there is one _ in the name
            underscored = '_'.join(target.split('_')[:-1]) if len(target.split('_')) > 2 else target.split('_')[0]
            try:
                # Since the number of core genes is the same for each reference strain, only need to determine it once
                if underscored == sorted(coregenomes)[0]:
                    coregenes.add(target)
                # If the percent identity is greater than the cutoff - adjust the cutoff to 90% for these analyses
                self.cutoff = 90
                if percentidentity >= self.cutoff:
                    # Update the dictionary with the target and the number of hits
                    resultdict[underscored] += 1
            except IndexError:
                pass
        # Sort the dictionary on the number of hits - best at the top
        topcore = sorted(resultdict.items(), key=operator.itemgetter(1), reverse=True)
        # Initialise a dictionary attribute to store results
        sample[self.analysistype].blastresults = dict()
        # If there are no results, populate negative results
        if not resultdict:
            sample[self.analysistype].blastresults = 'NA'
        # If results, add a string of the best number of hits, and a string of the total number of genes
        else:
            sample[self.analysistype].blastresults[topcore[0][0]] = (str(topcore[0][1]), str(len(coregenes)))

    def csvwriter(self):
        header = 'Strain,ClosestRef,GenesPresent/Total,\n'
        data = ''
        for sample in self.metadata:
            try:
                if sample[self.analysistype].blastresults != 'NA':
                    # Write the sample name, closest ref genome, and the number of genes found / total number of genes
                    closestref = sample[self.analysistype].blastresults.items()[0][0]
                    coregenes = sample[self.analysistype].blastresults.items()[0][1][0]
                    totalcore = sample[self.analysistype].blastresults.items()[0][1][1]
                    # Add the data to the object
                    sample[self.analysistype].targetspresent = coregenes
                    sample[self.analysistype].totaltargets = totalcore
                    row = '{},{},{}/{}'.format(sample.name, closestref, coregenes, totalcore)
                    # If the script is being run as part of the assembly pipeline, make a report for each sample
                    if self.pipeline:
                        # Open the report
                        with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                       self.analysistype), 'wb') as report:
                            # Write the row to the report
                            report.write(header)
                            report.write(row)
                    data += row
                    # Remove the messy blast results from the object
                    try:
                        delattr(sample[self.analysistype], "blastresults")
                    except KeyError:
                        pass
                else:
                    sample[self.analysistype].targetspresent = 'NA'
                    sample[self.analysistype].totaltargets = 'NA'
            except KeyError:
                sample[self.analysistype].targetspresent = 'NA'
                sample[self.analysistype].totaltargets = 'NA'
        with open('{}coregenome.csv'.format(self.reportpath), 'wb') as report:
            # Write the data to the report
            report.write(header)
            report.write(data)
