#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import *
from genesippr.genesippr import GeneSippr
__author__ = 'adamkoziol'


class Plasmids(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        header = 'Strain,Gene,PercentIdentity,Length,FoldCoverage\n'
        data = ''
        with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    multiple = False
                    for name, identity in sample[self.analysistype].results.items():
                        if not multiple:
                            data += '{},{},{},{}\n'.format(name, identity,
                                                           len(sample[self.analysistype].sequences[name]),
                                                           sample[self.analysistype].avgdepth[name])
                        else:
                            data += ',{},{},{},{}\n'.format(name, identity,
                                                            len(sample[self.analysistype].sequences[name]),
                                                            sample[self.analysistype].avgdepth[name])
                        multiple = True
                else:
                    data += '\n'
            report.write(header)
            report.write(data)


class Resistance(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create a set of all the gene names without alleles or accessions e.g. sul1_18_AY260546 becomes sul1
        genedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(self.targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            for line in notes:
                # Ignore comment lines - they will break the parsing
                if line.startswith('#'):
                    continue
                # Split the line on colons e.g. QnrB53:Quinolone resistance: has three variables after the split:
                # gene(QnrB53), resistance(Quinolone resistance), and _(\n)
                gene, resistance, _ = line.split(':')
                # Set up the resistance dictionary
                genedict[gene] = resistance
        # Find unique gene names with the highest percent identity
        for sample in self.runmetadata.samples:
            if sample[self.analysistype].results:
                # Initialise a dictionary to store the unique genes, and their percent identities
                sample[self.analysistype].uniquegenes = dict()
                for name, identity in sample[self.analysistype].results.items():
                    # Split the name of the gene from the string e.g. ARR-2_1_HQ141279 yields ARR-2
                    genename = name.split('_')[0]
                    # Set the best observed percent identity for each unique gene
                    try:
                        # Pull the previous best identity from the dictionary
                        bestidentity = sample[self.analysistype].uniquegenes[genename]
                        # If the current identity is better than the old identity, save it
                        if float(identity) > float(bestidentity):
                            sample[self.analysistype].uniquegenes[genename] = float(identity)
                    # Initialise the dictionary if necessary
                    except KeyError:
                        sample[self.analysistype].uniquegenes[genename] = float(identity)
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        header = 'Strain,Resistance,Gene,Allele,Accession,PercentIdentity,Length,FoldCoverage\n'
        data = str()
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    # If there are multiple results for a sample, don't write the sample name in each line of the report
                    multiple = False
                    for name, identity in sample[self.analysistype].results.items():
                        # Split the name on underscores: ARR-2_1_HQ141279; gene: ARR-2, allele: 1, accession: HQ141279
                        genename, allele, accession = name.split('_')
                        # Retrieve the best identity for each gene
                        percentid = sample[self.analysistype].uniquegenes[genename]
                        # If the percent identity of the current gene matches the best percent identity, add it to
                        # the report - there can be multiple occurrences of genes e.g.
                        # sul1,1,AY224185,100.00,840 and sul1,2,CP002151,100.00,927 are both included because they
                        # have the same 100% percent identity
                        if float(identity) == percentid:
                            # Treat the initial vs subsequent results for each sample slightly differently - instead
                            # of including the sample name, use an empty cell instead
                            if not multiple:
                                # Populate the results
                                data += '{},{},{},{},{},{},{}\n'.format(
                                    genedict[genename],
                                    genename,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                            else:
                                data += ',{},{},{},{},{},{},{}\n'.format(
                                    genedict[genename],
                                    genename,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                            multiple = True
                else:
                    data += '\n'
            # Write the strings to the file
            report.write(header)
            report.write(data)
