#!/usr/bin/env python
from threading import Thread
from accessoryFunctions import *
__author__ = 'adamkoziol'


class CHAS(object):

    def chas(self):
        """
        Call the necessary functions
        """
        # Determine primer binding
        self.primers()
        # Parse the ePCR results
        self.epcrparsethreads()
        # Create a report
        self.report()

    def primers(self):
        """Setup and create  threads for ePCR"""
        from glob import glob
        # Create the threads for the ePCR analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.epcr, args=())
                threads.setDaemon(True)
                threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                # Get the primers ready
                try:
                    sample[self.analysistype].primers = \
                        glob('{}geneseekr/{}/{}/primers/*.txt'.format(self.reffilepath, self.analysistype,
                                                                      sample.general.referencegenus))[0]
                    # Find the name of the probe file
                    sample[self.analysistype].probes = \
                        glob('{}geneseekr/{}/{}/probes/*.fa'.format(self.reffilepath, self.analysistype,
                                                                    sample.general.referencegenus))[0]
                    # Create the BLAST database of the probes (if necessary)
                    self.makeblastdb(sample[self.analysistype].probes)
                    # Initialise a list to store the names of the targets
                    sample[self.analysistype].targets = list()
                    # Open the primer file, and read the names of the targets into a list
                    with open(sample[self.analysistype].primers, 'rb') as primerfile:
                        for line in primerfile:
                            sample[self.analysistype].targets.append(line.split('\t')[0])
                # Organisms without primer/probe files will fail. Populate metadata with 'NA' values
                except IndexError:
                    sample[self.analysistype].primers = 'NA'
                    sample[self.analysistype].probes = 'NA'
                # Only try to process organisms with primer files
                if sample[self.analysistype].primers != 'NA':
                    # Make the output path
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                    make_path(sample[self.analysistype].reportdir)
                    # Set the base name of the output file
                    outfile = sample[self.analysistype].reportdir + sample.name
                    # Set the hashing and mapping commands
                    sample.commands.famap = 'famap -b {}.famap {}.fasta'.format(outfile, sample.general.filenoext)
                    sample.commands.fahash = 'fahash -b {}.hash {}.famap'.format(outfile, outfile)
                    # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
                    # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
                    # -m 10000 (Set variability for STS size for lookup),
                    # -n 1 (Set max allowed mismatches per primer for lookup)
                    # -g 0 (Set max allowed indels per primer for lookup),
                    # -G (Print alignments in comments), -o {output file}
                    sample.commands.epcr = 're-PCR -S {}.hash -r + -m 10000 -n 2 -g 0 -G -q -o {}.txt {}' \
                        .format(outfile, outfile, sample[self.analysistype].primers)
                    # Add the variables to the queue
                    self.epcrqueue.put((sample, outfile))
        self.epcrqueue.join()

    def epcr(self):
        from subprocess import call
        while True:
            sample, linkfile = self.epcrqueue.get()
            # Set the names of the output files
            sample[self.analysistype].famap = '{}.famap'.format(linkfile)
            sample[self.analysistype].hash = '{}.hash'.format(linkfile)
            sample[self.analysistype].output = '{}.txt'.format(linkfile)
            # Initialise a list to store the results
            sample[self.analysistype].epcrresults = list()
            # If the files created by the results do not exist, run the necessary system calls
            if not os.path.isfile(sample[self.analysistype].famap):
                call(sample.commands.famap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            if not os.path.isfile(sample[self.analysistype].hash):
                call(sample.commands.fahash, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            if not os.path.isfile(sample[self.analysistype].output):
                call(sample.commands.epcr, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            # Read the results into a list
            with open(sample[self.analysistype].output, 'rb') as results:
                for line in results:
                    sample[self.analysistype].epcrresults.append(line.strip())
            self.epcrqueue.task_done()

    def epcrparsethreads(self):
        """
        Parse the ePCR results, and run BLAST on the parsed results
        """
        from Bio import SeqIO
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.epcrparse, args=())
                threads.setDaemon(True)
                threads.start()
        for sample in self.metadata:
            if sample[self.analysistype].primers != 'NA':
                # Initialise a dictionary to store the SeqIO records of each assembly
                record = dict()
                # Initialise dictionaries to store results in the object
                sample[self.analysistype].blastresults = dict()
                sample[self.analysistype].rawblastresults = dict()
                # Load the records from the assembly into the dictionary
                for rec in SeqIO.parse(sample.general.bestassemblyfile, 'fasta'):
                    record[rec.id] = str(rec.seq)
                # Iterate through the ePCR results
                for line in sample[self.analysistype].epcrresults:
                    # The data of interest is in the lines that do not start with a #
                    # TLH	2016-SEQ-0359_4_length_321195_cov_28.6354_ID_3773	+	227879	228086	0	0	208/1000-1000
                    if not line.startswith('#'):
                        # Add the variables to the queue
                        self.epcrparsequeue.put((sample, record, line))
        self.epcrparsequeue.join()

    def epcrparse(self):
        """
        Run BLAST, and record results to the object
        """
        from Bio.Blast.Applications import NcbiblastnCommandline
        while True:
            sample, record, line = self.epcrparsequeue.get()
            # Split the data on tabs
            gene, chromosome, strand, start, end, m_match, gaps, act_len_exp_len = line.split('\t')
            # Extract the gene sequence from the contigs
            # The record dictionary has the contig name, and the sequence. Splice out the data using the start and
            # end coordinates specified by ePCR
            genesequence = record[chromosome][int(start) - 1:int(end)]
            # Set up BLASTn using blastn-short, as the probe sequences tend to be very short
            blastn = NcbiblastnCommandline(db=sample[self.analysistype].probes.split('.')[0],
                                           num_threads=12,
                                           task='blastn-short',
                                           num_alignments=1,
                                           outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                  "evalue bitscore slen length'")
            # Run the BLASTn, with the gene sequence as stdin
            out, err = blastn(stdin=genesequence)
            # Split the output string on tabs
            results = out.rstrip().split('\t')
            # Populate the raw blast results
            sample[self.analysistype].rawblastresults[gene] = results
            # Create named variables from the list
            positives = float(results[2])
            mismatches = float(results[3])
            gaps = float(results[4])
            subjectlength = float(results[7])
            # Calculate the percent identity
            percentidentity = float('{:0.2f}'.format((positives - gaps) / subjectlength * 100))
            # Create a dictionary with the desired values to store in the metadata object
            resultdict = {
                'matches': positives,
                'mismatches': mismatches,
                'gaps': gaps,
                'subject_length': subjectlength,
                'percent_identity': percentidentity,
                'match_length': results[8].split('\n')[0]
            }
            # Populate the metadata object with the dictionary
            sample[self.analysistype].blastresults[gene] = resultdict
            self.epcrparsequeue.task_done()

    def makeblastdb(self, fastapath):
        """
        Makes blast database files from targets as necessary
        """
        import subprocess
        import shlex
        # remove the path and the file extension for easier future globbing
        db = fastapath.split('.')[0]
        nhr = '{}.nhr'.format(db)  # add nhr for searching
        if not os.path.isfile(str(nhr)):  # if check for already existing dbs
            # Create the databases
            subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                        .format(fastapath, db)), stdout=self.fnull, stderr=self.fnull)
        dotter()

    def report(self):
        """
        Create reports of the findings
        """
        # Initialise a variable to store the results
        data = ''
        for sample in self.metadata:
            if sample[self.analysistype].primers != 'NA':
                # Set the name of the strain-specific report
                sample[self.analysistype].report = '{}/{}_{}.csv'.format(sample[self.analysistype].reportdir,
                                                                         sample.name, self.analysistype)
                # Populate the strain-specific string with header, and strain name
                strainspecific = 'Strain,{},\n{},'.format(','.join(sorted(sample[self.analysistype].targets)),
                                                          sample.name)
                # Iterate through all the genes in the organism-specific analysis
                for gene in sorted(sample[self.analysistype].targets):
                    try:
                        # Extract the percent identity
                        percentidentity = sample[self.analysistype].blastresults[gene]['percent_identity']
                        # If the % identity is greater than the cutoff of 50%, the gene is considered to be present
                        if percentidentity > 50:
                            strainspecific += '{},'.format(percentidentity)
                        else:
                            strainspecific += '-,'
                    # If there are no BLAST results, then the gene is absent
                    except KeyError:
                        strainspecific += '-,'
                strainspecific += '\n'
                # Open and write the data to the strain-specific report
                with open(sample[self.analysistype].report, 'w') as specificreport:
                    specificreport.write(strainspecific)
                # Add all the data from each strain to the cumulative data string
                data += strainspecific
        # Open and write the cumulative data to the cumulative report
        with open('{}/{}.csv'.format(self.reportdir, self.analysistype), 'w') as report:
            report.write(data)

    def __init__(self, inputobject, analysistype):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.path = inputobject.path
        self.start = inputobject.start
        self.reffilepath = inputobject.referencefilepath
        self.threads = inputobject.threads
        self.reportdir = '{}/'.format(inputobject.reportdir)
        self.analysistype = analysistype
        self.epcrqueue = Queue(maxsize=self.threads)
        self.epcrparsequeue = Queue(maxsize=self.threads)
        self.fnull = open(os.path.devnull, 'wb')
        # Run the analyses
        self.chas()
