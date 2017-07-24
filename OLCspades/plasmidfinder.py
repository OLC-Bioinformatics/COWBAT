#!/usr/bin/env python
from prophages import *
__author__ = 'adamkoziol'


class PlasmidFinder(Prophages):
    def reporter(self):
        combinedrow = ''
        for sample in self.metadata:
            # Populate the header with the appropriate data
            row = 'Contig,PlasmidAccession,PercentIdentity\n'
            try:
                for result in sample[self.analysistype].blastresults.items():
                    contig = list(result)[0]
                    query = list(result)[1].keys()[0]
                    percentid = list(result)[1].values()[0]
                    # Add the data to the row
                    row += '{},{},{}\n'.format(contig, query, percentid)
            except (KeyError, AttributeError):
                row = ''
            combinedrow += row
            # If the script is being run as part of the assembly pipeline, make a report for each sample
            if self.pipeline:
                if sample.general.bestassemblyfile != 'NA':
                    # Open the report
                    with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                   self.analysistype), 'w') as report:
                        # Write the row to the report
                        report.write(row)
            # Remove the messy blast results from the object
            try:
                delattr(sample[self.analysistype], "blastresults")
            except KeyError:
                pass
        # Create the report containing all the data from all samples
        with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'w') \
                as combinedreport:
            combinedreport.write(combinedrow)
