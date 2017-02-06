#!/usr/bin/env python
import shlex
import subprocess
import time
from collections import defaultdict
from csv import DictReader
from glob import glob
from threading import Thread
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from accessoryFunctions import *

__author__ = 'mike knowles, adamkoziol'


class GeneSeekr(object):

    def geneseekr(self):
        # Make blast databases (if necessary)
        printtime('Creating {} blast databases as required'.format(self.analysistype), self.start)
        self.makedbthreads()
        # Run the blast analyses
        printtime('Running {} blast analyses'.format(self.analysistype), self.start)
        self.blastnthreads()
        globalcounter()
        self.reporter()
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata:
            delattr(sample[self.analysistype], "targetnames")
            delattr(sample[self.analysistype], "targets")
        printtime('{} analyses complete'.format(self.analysistype), self.start)

    def makedbthreads(self):
        """
        Setup and create threads for class
        """
        # Find all the target folders in the analysis and add them to the targetfolders set
        for sample in self.metadata:
            if sample[self.analysistype].combinedtargets != 'NA':
                self.targetfolders.add(sample[self.analysistype].targetpath)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.targetfolders)):
            # Send the threads to makeblastdb
            threads = Thread(target=self.makeblastdb, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Make blast databases for MLST files (if necessary)
        for targetdir in self.targetfolders:
            # List comprehension to remove any previously created database files from list
            self.targetfiles = glob('{}/*.tfa'.format(targetdir))
            for targetfile in self.targetfiles:
                # Read the sequences from the target file to a dictionary
                self.records[targetfile] = SeqIO.to_dict(SeqIO.parse(targetfile, 'fasta'))
                # Add the fasta file to the queue
                self.dqueue.put(targetfile)
        self.dqueue.join()  # wait on the dqueue until everything has been processed

    def makeblastdb(self):
        """Makes blast database files from targets as necessary"""
        while True:  # while daemon
            fastapath = self.dqueue.get()  # grabs fastapath from dqueue
            # remove the path and the file extension for easier future globbing
            db = fastapath.split('.')[0]
            nhr = '{}.nhr'.format(db)  # add nhr for searching
            fnull = open(os.devnull, 'w')  # define /dev/null
            if not os.path.isfile(str(nhr)):  # if check for already existing dbs
                # Create the databases
                # TODO use MakeBLASTdb class
                subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                            .format(fastapath, db)), stdout=fnull, stderr=fnull)
            dotter()
            self.dqueue.task_done()  # signals to dqueue job is done

    def blastnthreads(self):
        """Setup and create  threads for blastn and xml path"""
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                for i in range(len(sample[self.analysistype].combinedtargets)):
                    threads = Thread(target=self.runblast, args=())
                    threads.setDaemon(True)
                    threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample[self.analysistype].combinedtargets != 'NA':
                # Add each fasta file combination to the threads
                self.blastqueue.put((sample.general.bestassemblyfile, sample[self.analysistype].combinedtargets,
                                     sample))
        # Join the threads
        self.blastqueue.join()

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            make_path(sample[self.analysistype].reportdir)
            try:
                report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
                size = os.path.getsize(report)
                # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
                # remove the report, and run the blast analyses again
                if size == 0:
                    os.remove(report)
                    report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                              time.strftime("%Y.%m.%d.%H.%M.%S"))
            except IndexError:
                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
            db = target.split('.')[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly, db=db, evalue='1E-5', num_alignments=1000000,
                                           num_threads=12,
                                           # outfmt="'6 qseqid sseqid positive mismatch gaps "
                                           #        "evalue bitscore slen length'",
                                           outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                  "evalue bitscore slen length qstart qend qseq sstart send'",
                                           out=report)
            # Save the blast command in the metadata
            sample[self.analysistype].blastcommand = str(blastn)
            # Only run blast if the report doesn't exist
            if not os.path.isfile(report):
                try:
                    blastn()
                except:
                    self.blastqueue.task_done()
                    self.blastqueue.join()
                    try:
                        os.remove(report)
                    except IOError:
                        pass
                    raise
            # Run the blast parsing module
            self.blastparser(report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        """
        Parse the blast results, and store necessary data in dictionaries in sample object
        :param report: Name of the blast output report being parsed
        :param sample: sample object
        """
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        resultdict = dict()
        # Initialise a dictionary to store all the target sequences
        sample[self.analysistype].targetsequence = dict()
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            target = row['subject_id']
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Update the dictionary with the target and percent identity
                resultdict.update({target: percentidentity})
                # Determine if the orientation of the sequence is reversed compared to the reference
                if int(row['subject_end']) < int(row['subject_start']):
                    # Create a sequence object using Biopython
                    seq = Seq(row['query_sequence'], IUPAC.unambiguous_dna)
                    # Calculate the reverse complement of the sequence
                    querysequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    querysequence = row['query_sequence']
                # Add the sequence in the correct orientation to the sample
                sample[self.analysistype].targetsequence[target] = querysequence
            # Add the percent identity to the object
            sample[self.analysistype].blastresults = resultdict
        # Populate missing results with 'NA' values
        if not resultdict:
            sample[self.analysistype].blastresults = 'NA'

    def reporter(self):
        """
        Creates .xlsx reports using xlsxwriter
        """
        import xlsxwriter
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook('{}/{}.xlsx'.format(self.reportpath, self.analysistype))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        for sample in self.metadata:
            # Reset the column to zero
            col = 0
            # Initialise a list to store all the data for each strain
            data = list()
            # Initialise a list of all the headers with 'Strain'
            headers = ['Strain']
            if sample[self.analysistype].targetnames != 'NA':
                # Append the sample name to the data list only if the script could find targets
                data.append(sample.name)
                if sample[self.analysistype].blastresults != 'NA':
                    for target in sorted(sample[self.analysistype].targetnames):
                        # Add the name of the gene to the header
                        headers.append(target)
                        if sample[self.analysistype].blastresults != 'NA':
                            try:
                                # Append the percent identity to the data list
                                data.append(str(sample[self.analysistype].blastresults[target]))
                                # Only if the alignment option is selected, for inexact results, add alignments
                                if self.align and sample[self.analysistype].blastresults[target] != 100.00:
                                    # Align the protein (and nucleotide) sequences to the reference
                                    self.alignprotein(sample, target)
                                    # Add the appropriate headers
                                    headers.extend(['{}_aa_Alignment'.format(target),
                                                    '{}_aa_SNP_location'.format(target),
                                                    '{}_nt_Alignment'.format(target),
                                                    '{}_nt_SNP_location'.format(target)
                                                    ])
                                    # Add the alignment, and the location of mismatches for both nucleotide and amino
                                    # acid sequences
                                    data.extend([sample[self.analysistype].aaalign[target],
                                                 sample[self.analysistype].aaindex[target],
                                                 sample[self.analysistype].ntalign[target],
                                                 sample[self.analysistype].ntindex[target],
                                                 ])
                            # If there are no blast results for the target, add a '-'
                            except (KeyError, TypeError):
                                data.append('-')
                        # If there are no blast results at all, add a '-'
                        else:
                            data.append('-')
            # Write the header to the spreadsheet
            for header in headers:
                worksheet.write(row, col, header, bold)
                # Set the column width based on the longest header
                try:
                    columnwidth[col] = len(header)if len(header) > columnwidth[col] else columnwidth[col]
                except KeyError:
                    columnwidth[col] = len(header)
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
            # Increment the row and reset the column to zero in preparation of writing results
            row += 1
            col = 0
            # List of the number of lines for each result
            totallines = list()
            # Write out the data to the spreadsheet
            for results in data:
                worksheet.write(row, col, results, courier)
                try:
                    # Counting the length of multi-line strings yields columns that are far too wide, only count
                    # the length of the string up to the first line break
                    alignmentcorrect = len(results.split('\n')[0])
                    # Count the number of lines for the data
                    lines = results.count('\n')
                    # Add the number of lines to the list
                    totallines.append(lines)
                # If there are no newline characters, set the width to the length of the string
                except AttributeError:
                    alignmentcorrect = len(results)
                # Increase the width of the current column, if necessary
                columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else columnwidth[col]
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
            # Set the width of the row to be the number of lines (number of newline characters) * 11
            worksheet.set_row(row, max(totallines) * 12)
            # Increase the row counter for the next strain's data
            row += 1
        # Close the workbook
        workbook.close()

    def alignprotein(self, sample, target):
        """
        Create alignments of the sample nucleotide and amino acid sequences to the reference sequences
        """
        import re
        # Initialise dictionaries
        sample[self.analysistype].dnaseq = dict()
        sample[self.analysistype].protseq = dict()
        sample[self.analysistype].ntindex = dict()
        sample[self.analysistype].aaindex = dict()
        sample[self.analysistype].ntalign = dict()
        sample[self.analysistype].aaalign = dict()
        # In order to properly translate the nucleotide sequence, BioPython requests that the sequence is a multiple of
        # three - not partial codons. Trim the sequence accordingly
        remainder = 0 - len(sample[self.analysistype].targetsequence[target]) % 3
        seq = sample[self.analysistype].targetsequence[target] if remainder == 0 \
            else sample[self.analysistype].targetsequence[target][:remainder]
        # Set the DNA and protein sequences of the target in the sample
        sample[self.analysistype].dnaseq[target] = Seq(seq, IUPAC.unambiguous_dna)
        # Translate the nucleotide sequence
        sample[self.analysistype].protseq[target] = str(sample[self.analysistype].dnaseq[target].translate())
        for targetfile in self.targetfiles:
            # Trim the reference sequence to multiples of three
            refremainder = 0 - len(self.records[targetfile][target].seq) % 3
            refseq = str(self.records[targetfile][target].seq) if refremainder % 3 == 0 \
                else str(self.records[targetfile][target].seq)[:refremainder]
            # Translate the nucleotide sequence of the reference sequence
            refdna = Seq(refseq, IUPAC.unambiguous_dna)
            refprot = str(refdna.translate())
            # Align the nucleotide sequence of the reference to the sample. If the corresponding bases match, add
            # a |, otherwise a space
            ntalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                      zip(sample[self.analysistype].dnaseq[target], refdna)))
            # Create the nucleotide alignment: the sample sequence, the (mis)matches, and the reference sequence
            sample[self.analysistype].ntalign[target] = self.interleaveblastresults(
                sample[self.analysistype].dnaseq[target],
                refdna)
            # Regex to determine location of mismatches in the sequences
            sample[self.analysistype].ntindex[target] = ';'.join(
                (str(snp.start()) for snp in re.finditer(' ', ntalignment)))
            # Perform the same steps, except for the amino acid sequence
            aaalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                      zip(sample[self.analysistype].protseq[target], refprot)))
            sample[self.analysistype].aaalign[target] = self.interleaveblastresults(
                sample[self.analysistype].protseq[target],
                refprot)
            sample[self.analysistype].aaindex[target] = ';'.join(
                (str(snp.start()) for snp in re.finditer(' ', aaalignment)))

    @staticmethod
    def interleaveblastresults(query, subject):
        """
        Creates an interleaved string that resembles BLAST sequence comparisons
        :param query: Query sequence
        :param subject: Subject sequence
        :return: Properly formatted BLAST-like sequence comparison
        """
        # Initialise strings to hold the matches, and the final BLAST-formatted string
        matchstring = ''
        blaststring = ''
        # Iterate through the query
        for i, bp in enumerate(query):
            # If the current base in the query is identical to the corresponding base in the reference, append a '|'
            # to the match string, otherwise, append a ' '
            if bp == subject[i]:
                matchstring += '|'
            else:
                matchstring += ' '
        # Set a variable to store the progress through the sequence
        prev = 0
        # Iterate through the query, from start to finish in steps of 60 bp
        for j in range(0, len(query), 60):
            # BLAST results string. The components are: current position (padded to four characters), 'OLC', query
            # sequence, \n, matches, \n, 'ref', subject sequence. Repeated until all the sequence data are present.
            """
            0000 OLC ATGAAGAAGATATTTGTAGCGGCTTTATTTGCTTTTGTTTCTGTTAATGCAATGGCAGCT
                     ||||||||||| ||| | |||| ||||||||| || ||||||||||||||||||||||||
                 ref ATGAAGAAGATGTTTATGGCGGTTTTATTTGCATTAGTTTCTGTTAATGCAATGGCAGCT
            0060 OLC GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
                     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                 ref GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
            """
            blaststring += '{} OLC {}\n         {}\n     ref {}\n' \
                .format('{:04d}'.format(j), query[prev:j + 60], matchstring[prev:j + 60], subject[prev:j + 60])
            # Update the progress variable
            prev = j + 60
        # Return the properly formatted string
        return blaststring

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = inputobject.start
        self.analysistype = inputobject.analysistype
        self.reportpath = inputobject.reportdir
        self.targetfolders = set()
        self.targetfiles = list()
        self.records = dict()
        self.pipeline = inputobject.pipeline
        self.referencefilepath = inputobject.referencefilepath
        self.cpus = inputobject.threads
        self.align = inputobject.align
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence', 'subject_start', 'subject_end']
        self.plusdict = defaultdict(make_dict)
        self.dqueue = Queue(maxsize=self.cpus)
        self.blastqueue = Queue(maxsize=self.cpus)
        self.geneseekr()


def sequencenames(contigsfile):
    """
    Takes a multifasta file and returns a list of sequence names
    :param contigsfile: multifasta of all sequences
    :return: list of all sequence names
    """
    sequences = list()
    for record in SeqIO.parse(open(contigsfile, "rU"), "fasta"):
        sequences.append(record.id)
    return sequences

if __name__ == '__main__':

    class Parser(object):

        def strainer(self):
            from accessoryFunctions import GenObject, MetadataObject
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
                # Create the .mlst attribute
                setattr(metadata, self.analysistype, GenObject())
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = sample
                metadata[self.analysistype].targets = self.targets
                metadata[self.analysistype].combinedtargets = self.combinedtargets
                metadata[self.analysistype].targetpath = self.targetpath
                # metadata[self.analysistype].targetnames = [os.path.split(x)[1].split('.')[0] for x in self.targets]
                metadata[self.analysistype].targetnames = sequencenames(self.combinedtargets)
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
            parser.add_argument('-c', '--cutoff',
                                type=int,
                                default=70, help='Threshold for maximum unique bacteria for a single antibiotic')
            parser.add_argument('-n', '--numthreads',
                                type=int,
                                default=24,
                                help='Specify number of threads')
            parser.add_argument('-a', '--align',
                                action='store_true',
                                help='Optionally output alignments of genes with less than 100% identity to reference '
                                     'genes. This alignment will use amino acid sequences for both query and reference')
            args = parser.parse_args()
            self.sequencepath = os.path.join(args.sequencepath, '')
            assert os.path.isdir(self.sequencepath), 'Cannot locate sequence path as specified: {}'\
                .format(self.sequencepath)
            self.targetpath = os.path.join(args.targetpath, '')
            assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}'\
                .format(self.targetpath)
            self.reportpath = os.path.join(args.reportpath, '')
            make_path(self.reportpath)
            assert os.path.isdir(self.reportpath), 'Cannot locate report path as specified: {}'\
                .format(self.reportpath)
            self.cutoff = args.cutoff
            self.threads = args.numthreads
            self.align = args.align
            self.strains = list()
            self.targets = list()
            self.combinedtargets = str()
            self.samples = list()
            self.analysistype = 'geneseekr'
            self.start = time.time()
            self.strainer()

    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.start = self.runmetadata.start
            self.analysistype = self.runmetadata.analysistype
            self.cutoff = self.runmetadata.cutoff
            self.threads = int(self.runmetadata.threads)
            self.reportdir = self.runmetadata.reportpath
            self.pipeline = False
            self.referencefilepath = str()
            self.align = self.runmetadata.align
            # Run the analyses
            GeneSeekr(self)

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                if self.genusspecific:
                    # Allow Shigella to use the same targets as Escherichia
                    genus = sample.general.referencegenus if sample.general.referencegenus != 'Shigella' \
                        else 'Escherichia'
                    targetpath = '{}{}/{}'.format(self.referencefilepath, self.analysistype, genus)
                else:
                    targetpath = '{}{}/'.format(self.referencefilepath, self.analysistype)
                targets = glob('{}/*.fa*'.format(targetpath))
                targetcheck = glob('{}/*.*fa*'.format(targetpath))
                if targetcheck:
                    try:
                        combinedtargets = glob('{}/*.tfa'.format(targetpath))[0]
                    except IndexError:
                        combinetargets(targets, targetpath)
                        combinedtargets = glob('{}/*.tfa'.format(targetpath))[0]
                    sample[self.analysistype].targets = targets
                    sample[self.analysistype].combinedtargets = combinedtargets
                    sample[self.analysistype].targetpath = targetpath
                    sample[self.analysistype].targetnames = sequencenames(combinedtargets)
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                else:
                    # Set the metadata file appropriately
                    sample[self.analysistype].targets = 'NA'
                    sample[self.analysistype].combinedtargets = 'NA'
                    sample[self.analysistype].targetpath = 'NA'
                    sample[self.analysistype].targetnames = 'NA'
                    sample[self.analysistype].reportdir = 'NA'
                # Special typing for Vibrio involves in silico qPCR primer/probe binding
                if sample.general.referencegenus == 'Vibrio':
                    self.chas.append(sample)
            else:
                # Set the metadata file appropriately
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].targets = 'NA'
                sample[self.analysistype].combinedtargets = 'NA'
                sample[self.analysistype].targetpath = 'NA'
                sample[self.analysistype].targetnames = 'NA'
                sample[self.analysistype].reportdir = 'NA'
        if self.chas:
            import CHAS
            CHAS.CHAS(self, 'chas')

    def __init__(self, inputobject, analysistype, genusspecific, cutoff):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.referencefilepath = inputobject.reffilepath
        self.threads = inputobject.cpus
        self.reportdir = '{}/'.format(inputobject.reportpath)
        self.cutoff = cutoff
        self.pipeline = True
        self.genusspecific = genusspecific
        self.chas = list()
        self.align = False
        # Get the alleles and profile into the metadata
        self.strainer()
