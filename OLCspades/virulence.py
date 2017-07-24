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
            try:
                # Virulence schemes differ depending on the organism
                if sample.general.referencegenus == 'Escherichia':
                    row += 'Contig,Gene,Percentidentity\n'
                    geneset = sorted(set([gene[0].split(':')[0] for gene in
                                          sample[self.analysistype].blastresults.items()]))
                    for gene in geneset:
                        besthits = sorted([item for item in sample[self.analysistype].blastresults.items()
                                           if list(item)[0].startswith(gene)], key=lambda x: x[1], reverse=True)
                        if 'stx' in gene:
                            for hit in besthits:
                                if hit[1] > 98 and '{}{}'.format(gene, hit[0].split(',')[0].split(':')[-1]) not in row:
                                    # print hit, hit[0].split(',')[0].split(':')[-1]
                                    if len(hit[0].split(',')[0].split(':')[-1]) == 1:
                                        row += '{},{}{},{}\n'\
                                            .format(hit[0].split(',')[1], gene, hit[0].split(',')[0]
                                                    .split(':')[-1], hit[1])
                                    elif gene not in row:
                                        row += '{},{},{}\n'.format(hit[0].split(',')[1], gene, hit[1])
                            sample.general.stx = True
                            # print besthits[0][1]
                        else:
                            row += '{},{},{}\n'.format(besthits[0][0].split(',')[1], gene, besthits[0][1])
                            try:
                                if not sample.general.stx:
                                    sample.general.stx = False
                            except KeyError:
                                sample.general.stx = False
                else:
                    sample.general.stx = False
                combinedrow += row
                # If the script is being run as part of the assembly pipeline, make a report for each sample
                if self.pipeline:
                    if sample.general.bestassemblyfile != 'NA' and sample[self.analysistype].reportdir != 'NA':
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

# Run stand alone
if __name__ == '__main__':
    class Parser(object):

        def strainer(self):
            from .accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = sorted(glob('{}*.fa*'.format(self.sequencepath)))
            self.targets = sorted(glob('{}*.fa*'.format(self.targetpath)))
            try:
                self.combinedtargets = glob('{}/*.tfa'.format(self.targetpath))[0]
            except IndexError:
                combinetargets(self.targets, self.targetpath)
                self.combinedtargets = glob('{}/*.tfa'.format(self.targetpath))[0]
            # Populate the metadata object. This object will be populated to mirror the objects created in the
            # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
            # of a pipeline
            assert self.strains, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your sequence path is correct'.format(self.sequencepath)
            assert self.targets, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your target path is correct'.format(self.targetpath)
            for sample in self.strains:
                # Create the object
                metadata = MetadataObject()
                # Set the base file name of the sequence. Just remove the file extension
                filename = os.path.split(sample)[1].split('.')[0]
                # Set the .name attribute to be the file name
                metadata.name = filename
                # Create the .general attribute
                metadata.general = GenObject()
                metadata.general.referencegenus = 'Escherichia'
                # Create the .mlst attribute
                setattr(metadata, self.analysistype, GenObject())
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = sample
                metadata[self.analysistype].targets = self.targets
                metadata[self.analysistype].combinedtargets = self.combinedtargets
                metadata[self.analysistype].targetpath = self.targetpath
                metadata[self.analysistype].targetnames = [os.path.split(x)[1].split('.')[0] for x in self.targets]
                metadata[self.analysistype].reportdir = self.reportpath
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            parser = ArgumentParser(description='Use to find markers for any bacterial genome')
            parser.add_argument('--version',
                                action='version',
                                version='%(prog)s v0.5')
            parser.add_argument('-s', '--sequencepath',
                                required=True,
                                help='Specify input fasta folder')
            parser.add_argument('-t', '--targetpath',
                                required=True,
                                help='Specify folder of targets')
            parser.add_argument('-r', '--reportpath',
                                required=True,
                                help='Specify output folder for csv')
            # parser.add_argument('-a', '--anti', type=str, required=True, help='JSON file location')
            parser.add_argument('-c', '--cutoff',
                                type=int,
                                default=70, help='Threshold for maximum unique bacteria for a single antibiotic')
            parser.add_argument('-n', '--numthreads',
                                type=int,
                                default=24,
                                help='Specify number of threads')
            args = parser.parse_args()
            self.sequencepath = os.path.join(args.sequencepath, '')
            assert os.path.isdir(self.sequencepath), 'Cannot locate sequence path as specified: {}' \
                .format(self.sequencepath)
            self.targetpath = os.path.join(args.targetpath, '')
            assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}' \
                .format(self.targetpath)
            self.reportpath = os.path.join(args.reportpath, '')
            assert os.path.isdir(self.reportpath), 'Cannot locate report path as specified: {}' \
                .format(self.reportpath)
            self.cutoff = args.cutoff
            self.threads = args.numthreads
            #
            self.strains = []
            self.targets = []
            self.combinedtargets = ''
            self.samples = []
            self.analysistype = 'virulence'
            self.start = time.time()
            self.strainer()


    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.start = self.runmetadata.start
            self.analysistype = self.runmetadata.analysistype
            # self.alleles = self.runmetadata.alleles
            # self.profile = self.runmetadata.profile
            self.cutoff = self.runmetadata.cutoff
            self.threads = self.runmetadata.threads
            self.reportdir = self.runmetadata.reportpath
            self.pipeline = False
            self.referencefilepath = ''
            # Run the analyses
            Virulence(self)

    # Run the class
    MetadataInit()
