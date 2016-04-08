#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class Prophages(GeneSeekr):

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
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Update the dictionary with the target and percent identity
                resultdict.update({row['query_id']: {row['subject_id']: percentidentity}})
            sample[self.analysistype].blastresults = resultdict

    def csvwriter(self):
        from csv import DictReader
        # Set the required variables to load prophage data from a summary file
        targetpath = '{}{}/'.format(self.referencefilepath, self.analysistype)
        overview = glob('{}/*.txt'.format(targetpath))[0]
        fieldnames = ['id_prophage', 'file_name', 'host', 'host_id', 'number_of_prophages_in_host',
                      'start_position_of_prophage', 'end_position_of_prophage', 'length_of_prophage']
        combinedrow = ''
        for sample in self.metadata:
            try:
                # Populate the header with the appropriate data
                row = 'Contig,ProphageID,Host,PercentIdentity\n'
                for result in sample[self.analysistype].blastresults.items():
                    # Open the prophage file as a dictionary - I do this here, as if I open it earlier, it looks like the
                    # file remains partially-read through for the next iteration. Something like prophagedata.seek(0) would
                    # probably work, but Dictreader objects don't have a .seek attribute
                    prophagedata = DictReader(open(overview), fieldnames=fieldnames, dialect='excel-tab')
                    # Set variable names for the unsightly stored values
                    contig = list(result)[0]
                    query = list(result)[1].keys()[0]
                    percentid = list(result)[1].values()[0]
                    # Iterate through the phage data in the dictionary
                    for phage in prophagedata:
                        if phage['id_prophage'] == query:
                            # Add the data to the row
                            row += '{},{},{},{}\n'.format(contig, query, phage['host'], percentid)
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
            except (AttributeError, KeyError):
                pass
        # Create the report containing all the data from all samples
        with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                as combinedreport:
            combinedreport.write(combinedrow)
