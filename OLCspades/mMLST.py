#!/usr/bin/env python
import json
import operator
import re
import shlex
import shutil
import subprocess
import time
from Queue import Queue
from collections import defaultdict
from csv import DictReader
from glob import glob
from threading import Thread

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

import getmlst
from accessoryFunctions import *

__author__ = 'akoziol, mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
Revised with speed improvements
"""


class MLST(object):
    def mlst(self):
        # Get the MLST profiles into a dictionary for each sample
        printtime('Populating {} sequence profiles'.format(self.analysistype), self.start)
        self.profiler()
        globalcounter()
        # Make blast databases (if necessary)
        printtime('Creating {} blast databases as required'.format(self.analysistype), self.start)
        self.makedbthreads(self.allelefolders)
        # Run the blast analyses
        printtime('Running {} blast analyses'.format(self.analysistype), self.start)
        self.blastnprep()
        globalcounter()
        # Determine sequence types from the analyses
        printtime('Determining {} sequence types'.format(self.analysistype), self.start)
        self.sequencetyper()
        globalcounter()
        # Optionally dump :self.resultprofile to :self.reportpath
        if self.datadump:
            self.dumper()
            printtime('{} reference profile dump complete'.format(self.analysistype), self.start)
        # Optionally determine the closest reference genome from a pre-computed profile (this profile would have been
        # created using self.datadump
        if self.bestreferencegenome and self.analysistype.lower() == 'rmlst':
            printtime('Finding closest reference genomes'.format(self.analysistype), self.start)
            self.referencegenomefinder()
        # Create reports
        printtime('Creating {} reports'.format(self.analysistype), self.start)
        self.reporter()
        globalcounter()
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata:
            try:
                delattr(sample[self.analysistype], "allelenames")
                delattr(sample[self.analysistype], "alleles")
                delattr(sample[self.analysistype], "profiledata")
            except KeyError:
                pass
        printtime('{} analyses complete'.format(self.analysistype), self.start)

    def profiler(self):
        """Creates a dictionary from the profile scheme(s)"""
        # Initialise variables
        profiledata = defaultdict(make_dict)
        profileset = set()
        supplementalset = ''
        genedict = {}
        # Find all the unique profiles to use with a set
        for sample in self.metadata:
            if sample[self.analysistype].profile != 'NA':
                profileset.add(sample[self.analysistype].profile[0])
                if self.analysistype == 'rmlst':
                    supplementalset = sample[self.analysistype].supplementalprofile

        # Extract the profiles for each set
        for sequenceprofile in profileset:
            # Clear the list of genes
            genelist = []
            for sample in self.metadata:
                if sequenceprofile == sample[self.analysistype].profile[0]:
                    genelist = [os.path.split(x)[1].split('.')[0] for x in sample[self.analysistype].alleles]
            try:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile), dialect='excel-tab')
            # Revert to standard comma separated values
            except KeyError:
                # Open the sequence profile file as a dictionary
                profile = DictReader(open(sequenceprofile))
            # Iterate through the rows
            for row in profile:
                # Iterate through the genes
                for gene in genelist:
                    # Add the sequence profile, and type, the gene name and the allele number to the dictionary
                    try:
                        profiledata[sequenceprofile][row['ST']][gene] = row[gene]
                    except KeyError:
                        try:
                            profiledata[sequenceprofile][row['rST']][gene] = row[gene]
                        except KeyError:
                            raise

            # Load the supplemental profile definitions
            if self.analysistype == 'rmlst':
                supplementalprofile = DictReader(open(supplementalset), dialect='excel-tab')
                # Do the same with the supplemental profile
                for row in supplementalprofile:
                    # Iterate through the genes
                    for gene in genelist:
                        # Add the sequence profile, and type, the gene name and the allele number to the dictionary
                        profiledata[sequenceprofile][row['rST']][gene] = row[gene]
            # Add the gene list to a dictionary
            genedict[sequenceprofile] = sorted(genelist)
            # Add the profile data, and gene list to each sample
            for sample in self.metadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sequenceprofile == sample[self.analysistype].profile[0]:
                        # Populate the metadata with the profile data
                        sample[self.analysistype].profiledata = profiledata[sample[self.analysistype].profile[0]]
                        # Add the allele directory to a list of directories used in this analysis
                        self.allelefolders.add(sample[self.analysistype].alleledir)
                        dotter()

    def makedbthreads(self, folder):
        """
        Setup and create threads for class
        :param folder: folder with sequence files with which to create blast databases
        """
        # Create and start threads for each fasta file in the list
        for i in range(len(folder)):
            # Send the threads to makeblastdb
            threads = Thread(target=self.makeblastdb, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Make blast databases for MLST files (if necessary)
        for alleledir in folder:
            # List comprehension to remove any previously created database files from list
            allelefiles = glob('{}/*.fasta'.format(alleledir))
            # For each allele file
            for allelefile in allelefiles:
                # Add the fasta file to the queue
                self.dqueue.put(allelefile)
        self.dqueue.join()  # wait on the dqueue until everything has been processed

    def makeblastdb(self):
        """Makes blast database files from targets as necessary"""
        while True:  # while daemon
            fastapath = self.dqueue.get()  # grabs fastapath from dqueue
            # remove the path and the file extension for easier future globbing
            db = fastapath.split('.')[0]
            nhr = '{}.nhr'.format(db)  # add nhr for searching
            if not os.path.isfile(str(nhr)):  # if check for already existing dbs
                # Create the databases
                # TODO use MakeBLASTdb class
                subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                            .format(fastapath, db)), stdout=self.fnull, stderr=self.fnull)
            dotter()
            self.dqueue.task_done()  # signals to dqueue job is done

    def blastnprep(self):
        """Setup blastn analyses"""
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    for allele in sample[self.analysistype].combinedalleles:
                        # Add each fasta/allele file combination to the threads
                        self.runblast(sample.general.bestassemblyfile, allele, sample)

    def runblast(self, assembly, allele, sample):
        """
        Run the BLAST analyses
        :param assembly: assembly path/file
        :param allele: combined allele file
        :param sample: sample object
        :return:
        """
        genome = os.path.split(assembly)[1].split('.')[0]
        # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
        # Do not re-perform the BLAST search each time
        make_path(sample[self.analysistype].reportdir)
        try:
            report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
            size = os.path.getsize(report)
            if size == 0:
                os.remove(report)
                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
        except IndexError:
            report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                      time.strftime("%Y.%m.%d.%H.%M.%S"))
        db = allele.split('.')[0]
        # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        blastn = NcbiblastnCommandline(query=assembly, db=db, evalue='1E-20', num_alignments=1000000,
                                       num_threads=12,
                                       outfmt="'6 qseqid sseqid positive mismatch gaps "
                                              "evalue bitscore slen length'",
                                       out=report)
        # Save the blast command in the metadata
        sample[self.analysistype].blastcommand = str(blastn)
        if not os.path.isfile(report):
            # Run BLAST
            blastn()
        # Run the blast parsing module
        self.blastparser(report, sample)

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            if row['subject_length'] is not None:
                percentidentity = (float(row['positives']) - float(row['gaps'])) / float(row['subject_length']) * 100
                bitscore = float(row['bit_score'])
            # Find the allele number and the text before the number for different formats
            allelenumber, gene = allelesplitter(row['subject_id'])
            # If the percent identity is 100, and there are no mismatches, the allele is a perfect match
            if percentidentity == 100 and float(row['mismatches']) == 0:
                # If there are multiple best hits, then the .values() will be populated
                if self.plusdict[sample.name][gene].values():
                    # If the previous best hit have under 100% identity, or if the current bitscore is better
                    if self.plusdict[sample.name][gene].values()[0].keys()[0] < 100:
                        # Clear the previous match
                        self.plusdict[sample.name][gene].clear()
                        # Populate the new match
                        self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore
                    # If the bitscore is better (longer match) clear the previous result
                    # (not for rMLST analyses, which are allowed multiple allele matches)
                    else:
                        if bitscore > self.plusdict[sample.name][gene].values()[0].values() and \
                                        self.analysistype.lower() != 'rmlst':
                            # Clear the previous match
                            self.plusdict[sample.name][gene].clear()
                            # Populate the new match
                            self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore
                        else:
                            # Add the allele to the gene match
                            self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore

                # Populate the match
                else:
                    self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore
            # If the match is above the cutoff, but below 100%, add it to the dictionary
            elif percentidentity > self.cutoff:
                # If there are multiple best hits, then the .values() will be populated
                if self.plusdict[sample.name][gene].values():
                    if bitscore > self.plusdict[sample.name][gene].values()[0].values() and \
                                    self.plusdict[sample.name][gene].values()[0].keys()[0] < 100:
                        # elif percentidentity > self.cutoff and gene not in self.plusdict[sample.name]:
                        self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore
                else:
                    self.plusdict[sample.name][gene][allelenumber][percentidentity] = bitscore
        # Update alleles (if desired)
        if self.updateallele and self.analysistype != 'mlst':
            # Initialise variables
            previousallele, percidentity, bitscore, clearallele, geneofinterest = '', '', '', '', ''
            # Iterate through the all the genes in the analyses
            for gene in sample[self.analysistype].allelenames:
                # Find the allele and percent identity + bit score (percent score)
                for allele, percentscore in self.plusdict[sample.name][gene].items():
                    # Split the bit score from the percent identity
                    percentidentity = percentscore.items()[0][0]
                    # If the percent identity is less than 100%, run the allele updater method
                    if 90 < float(percentidentity) < 100:
                        geneofinterest, previousallele, percidentity, bitscore = \
                            self.alleleupdater(sample, gene, allele)
                self.plusdict[sample.name][geneofinterest].clear()
                self.plusdict[sample.name][geneofinterest][previousallele][percidentity] = bitscore
        # Populate empty results for genes without any matches
        for gene in sample[self.analysistype].allelenames:
            if gene not in self.plusdict[sample.name]:
                self.plusdict[sample.name][gene]['N'][0] = 0

    def alleleupdater(self, sample, gene, targetallele):
        """
        Updates file of alleles if the new allele passes length and identity checks
        :param sample: sample object
        :param gene: name of gene of interest
        :param targetallele: closest allele in database
        """
        from cStringIO import StringIO
        from Bio.Blast import NCBIXML
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        # Find the name of the allele file based on whether the gene name is found in the file name
        genefile = [x for x in sample[self.analysistype].alleles if gene in x][0]
        # Remove the extension of the gene file for use in makeblastdb and blastn
        genefilenoext = genefile.split(".")[0]
        # Create the blast databases every time the method is called
        subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                    .format(genefile, genefilenoext)), stdout=self.fnull, stderr=self.fnull)
        # Re-perform BLAST analyses, but using an XML output that stores the alignment
        blastn = NcbiblastnCommandline(query=sample.general.bestassemblyfile, db=genefilenoext, evalue=0.1, outfmt=5)
        # Note that there is no output file specified -  the search results are currently stored in stdout
        stdout, stderr = blastn()
        if stdout.find('Hsp') != -1:
            # Map the blast record to memory
            blast_handle = StringIO(stdout)
            # Open record from memory-mapped file
            records = NCBIXML.parse(blast_handle)
            # Iterate through the records in the blast results
            for record in records:  # This process is just to retrieve HSPs from xml files
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        # Extract the allele name from the blast result
                        allele = str(alignment.accession.split("_")[-1]) if "_" in alignment.accession else \
                            str(alignment.accession.split("-")[-1])
                        # The point of this blast is to extract the sequence of the allele
                        # If the allele name matches the target allele (the closest match e.g. BACT000063_20 -> 20)
                        if allele == str(targetallele):
                            # As there is some discrepancy with the capitalisation of terms, make sure it is consistent
                            analysistype = 'rMLST' if self.analysistype == 'rmlst' else 'MLST'
                            # Set the directory containing the profile and alleles
                            alleledir = self.referencefilepath + analysistype
                            # The name of the supplemental allele file (without and with the .fa extension)
                            allelefilenoext = '{}/OLC_{}_alleles'.format(alleledir, analysistype)
                            allelefile = allelefilenoext + '.fa'
                            # Create the file if it doesn't exist
                            open(allelefile, 'ab').close()
                            # Create a list of all the blast database files in the folder
                            dbfiles = glob('{}.n*'.format(allelefilenoext))
                            # Remove the database files
                            map(lambda y: os.remove(y), dbfiles)
                            # Create the necessary blast database files
                            subprocess.call(
                                shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                                            .format(allelefile, allelefilenoext)), stdout=self.fnull, stderr=self.fnull)

                            # Perform BLAST analysis using the supplemental allele file as the database
                            nestedblastn = NcbiblastnCommandline(db=allelefilenoext, evalue=0.1, outfmt=5)
                            # There is no output file specified; the search results are currently stored in stdout
                            # Additionally, the sequence from the previous BLAST query is used as stdin in this BLAST
                            import Bio.Application
                            try:
                                nestedstdout, nestedstderr = nestedblastn(stdin=hsp.query)
                            except Bio.Application.ApplicationError:
                                nestedstdout = ''
                                pass
                            # Search stdout for matches - if the term Hsp appears (the .find function will NOT
                            # return -1), a match has been found, and stdout is written to file
                            if nestedstdout.find('Hsp') != -1:
                                nested_blast_handle = StringIO(nestedstdout)
                                # Open record from memory-mapped file
                                nestedrecords = NCBIXML.parse(nested_blast_handle)
                                # Initialise variables
                                previousallele = ''
                                bitscore = ''
                                # Iterate through the records
                                for nestedrecord in nestedrecords:
                                    for nestedalignment in nestedrecord.alignments:
                                        for nestedhsp in nestedalignment.hsps:
                                            # Calculate the percent identity
                                            percentidentity = float("%.2f" % float(float(nestedhsp.identities) /
                                                                                   float(nestedalignment.length) * 100))
                                            # Only set the previous allele and the bitscore if the % identity is 100%
                                            if percentidentity == 100:
                                                # Get the allele number
                                                previousallele = nestedalignment.accession.split("_")[-1]
                                                bitscore = nestedhsp.score

                                # If this allele has already been found in the supplemental file, return these data
                                if previousallele:
                                    return gene, previousallele, 100.0, bitscore
                                # Otherwise, get the allele sequence into the supplemental file
                                else:
                                    # Initialise variables
                                    allelelist = []
                                    allelenumber = 1000000
                                    # Open the allele file to append
                                    with open(allelefile, 'ab+') as supplemental:
                                        for line in supplemental:
                                            # As there are multiple genes (e.g. BACT000001, BACT000063, etc.) check to
                                            # see if the gene of interest is in the line
                                            if gene in line:
                                                # Append the header to the list
                                                allelelist.append(line)
                                        # Find the last allele number from the header
                                        try:
                                            allelenumber = int(allelelist[-1].split("_")[1]) + 1
                                        # If there are no alleles for the gene of interest, then pass
                                        except IndexError:
                                            pass
                                        # Puts the HSP in the correct order -  hits to the negative strand will be
                                        # reversed compared to what we're looking for
                                        allelesequence = Seq(hsp.query, generic_dna)
                                        if hsp.sbjct_start < hsp.sbjct_end:
                                            end = hsp.sbjct_end
                                        else:
                                            end = hsp.sbjct_start
                                            allelesequence = allelesequence.reverse_complement()
                                        # Screen out hits that are shorter than the targets
                                        # Keeping this format though this statement could be re-written more efficiently
                                        if end < alignment.length:
                                            pass
                                        # The header will be >BACT00001_1000000
                                        definitionline = '{}_{}'.format(gene, allelenumber)
                                        # Create a sequence record using BioPython
                                        fasta = SeqRecord(allelesequence,
                                                          # Without this, the header will be improperly formatted
                                                          description='',
                                                          # Use >:definitionline as the header
                                                          id=definitionline)
                                        # Use the SeqIO module to properly format the new sequence record
                                        SeqIO.write(fasta, supplemental, "fasta")
                                        # Return the necessary values
                                        return gene, allelenumber, 100.0, hsp.score
                            # If there are hits in the supplemental allele file
                            else:
                                # Initialise variables
                                allelelist = []
                                allelenumber = 1000000
                                # Open the allele file to find the last allele associated with the gene of interest
                                with open(allelefile, 'ab+') as supplemental:
                                    for line in supplemental:
                                        if gene in line:
                                            allelelist.append(line)
                                    try:
                                        allelenumber = int(allelelist[-1].split("_")[1]) + 1
                                    except IndexError:
                                        pass
                                    # Puts the HSP in the correct order -  hits to the negative strand will be
                                    # reversed compared to what we're looking for
                                    allelesequence = Seq(hsp.query, generic_dna)
                                    if hsp.sbjct_start < hsp.sbjct_end:
                                        end = hsp.sbjct_end
                                    else:
                                        end = hsp.sbjct_start
                                        allelesequence = allelesequence.reverse_complement()
                                    # Screen out hits that are shorter than the targets
                                    # Keeping this format, though this if statement could be re-written more efficiently
                                    if end < alignment.length:
                                        pass
                                    definitionline = '{}_{}'.format(gene, allelenumber)
                                    # Create a sequence record using BioPython
                                    fasta = SeqRecord(allelesequence,
                                                      # If this is not added, the header will not be formatted properly
                                                      description='',
                                                      # Use >:definitionline as the header
                                                      id=definitionline)
                                    # Use the SeqIO module to properly format the new sequence record
                                    SeqIO.write(fasta, supplemental, "fasta")
                                    # Return the appropriate information
                                    return gene, allelenumber, 100.0, hsp.score

    def sequencetyper(self):
        """Determines the sequence type of each strain based on comparisons to sequence type profiles"""
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    # Initialise variables
                    header = 0
                    # Iterate through the genomes
                    # for sample in self.metadata:
                    genome = sample.name
                    # Initialise self.bestmatch[genome] with an int that will eventually be replaced by the # of matches
                    self.bestmatch[genome] = defaultdict(int)
                    if sample[self.analysistype].profile != 'NA':
                        # Create the profiledata variable to avoid writing self.profiledata[self.analysistype]
                        profiledata = sample[self.analysistype].profiledata
                        # For each gene in plusdict[genome]
                        for gene in sample[self.analysistype].allelenames:
                            # Clear the appropriate count and lists
                            multiallele = []
                            multipercent = []
                            # Go through the alleles in plusdict
                            for allele in self.plusdict[genome][gene]:
                                percentid = self.plusdict[genome][gene][allele].keys()[0]
                                # "N" alleles screw up the allele splitter function
                                if allele != "N":
                                    # Use the alleleSplitter function to get the allele number
                                    # allelenumber, alleleprenumber = allelesplitter(allele)
                                    # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                                    multiallele.append(int(allele))
                                    multipercent.append(percentid)
                                # If the allele is "N"
                                else:
                                    # Append "N" and a percent identity of 0
                                    multiallele.append("N")
                                    multipercent.append(0)
                                if not multiallele:
                                    multiallele.append("N")
                                    multipercent.append(0)
                            if self.analysistype == 'rmlst':
                                # For whatever reason, the rMLST profile scheme treat multiple allele hits as 'N's.
                                multiallele = multiallele if len(multiallele) == 1 else ['N']
                                if multipercent:
                                    multipercent = multipercent if len(multiallele) == 1 else [0, 0]
                                else:
                                    multipercent = [0]
                            # Populate self.bestdict with genome, gene, alleles joined with a space (this was made like
                            # this because allele is a list generated by the .iteritems() above
                            self.bestdict[genome][gene][" ".join(str(allele)
                                                                 for allele in sorted(multiallele))] = multipercent[0]
                            # Find the profile with the most alleles in common with the query genome
                            for sequencetype in profiledata:
                                # The number of genes in the analysis
                                header = len(profiledata[sequencetype])
                                # refallele is the allele number of the sequence type
                                refallele = profiledata[sequencetype][gene]
                                # If there are multiple allele matches for a gene in the reference profile e.g. 10 692
                                if len(refallele.split(" ")) > 1:
                                    # Map the split (on a space) alleles as integers - if they are treated as integers,
                                    # the alleles will sort properly
                                    intrefallele = map(int, refallele.split(" "))
                                    # Create a string of the joined, sorted alleles
                                    sortedrefallele = " ".join(str(allele) for allele in sorted(intrefallele))
                                else:
                                    # Use the reference allele as the sortedRefAllele
                                    sortedrefallele = refallele
                                for allele, percentid in self.bestdict[genome][gene].iteritems():
                                    # If the allele in the query genome matches the allele in the reference profile, add
                                    # the result to the bestmatch dictionary. Genes with multiple alleles were sorted
                                    # the same, strings with multiple alleles will match: 10 692 will never be 692 10
                                    if allele == sortedrefallele and float(percentid) == 100.00:
                                        # Increment the number of matches to each profile
                                        self.bestmatch[genome][sequencetype] += 1
                                    # Special handling of BACT000060 and BACT000065 genes. When the reference profile
                                    # has an allele of 'N', and the query allele doesn't, set the allele to 'N', and
                                    # count it as a match
                                    elif gene == 'BACT000060' or gene == 'BACT000065':
                                        if sortedrefallele == 'N' and allele != 'N':
                                            # Increment the number of matches to each profile
                                            self.bestmatch[genome][sequencetype] += 1
                                    elif allele == sortedrefallele and sortedrefallele == 'N':
                                        # Increment the number of matches to each profile
                                        self.bestmatch[genome][sequencetype] += 1
                        # Get the best number of matches
                        # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
                        try:
                            sortedmatches = sorted(self.bestmatch[genome].items(), key=operator.itemgetter(1),
                                                   reverse=True)[0][1]
                        # If there are no matches, set :sortedmatches to zero
                        except IndexError:
                            sortedmatches = 0
                        # Otherwise, the query profile matches the reference profile
                        if int(sortedmatches) == header:
                            # Iterate through best match
                            for sequencetype, matches in self.bestmatch[genome].iteritems():
                                if matches == sortedmatches:
                                    for gene in profiledata[sequencetype]:
                                        # Populate resultProfile with the genome, best match to profile, # of matches
                                        # to the profile, gene, query allele(s), reference allele(s), and % identity
                                        self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                            self.bestdict[genome][gene]
                                                .keys()[0]] = str(self.bestdict[genome][gene].values()[0])
                                    sample[self.analysistype].sequencetype = sequencetype
                                    sample[self.analysistype].matchestosequencetype = matches
                        # If there are fewer matches than the total number of genes in the typing scheme
                        elif 0 < int(sortedmatches) < header:
                            mismatches = []
                            # Iterate through the sequence types and the number of matches in bestDict for each genome
                            for sequencetype, matches in self.bestmatch[genome].iteritems():
                                # If the number of matches for a profile matches the best number of matches
                                if matches == sortedmatches:
                                    # Iterate through the gene in the analysis
                                    for gene in profiledata[sequencetype]:
                                        # Get the reference allele as above
                                        refallele = profiledata[sequencetype][gene]
                                        # As above get the reference allele split and ordered as necessary
                                        if len(refallele.split(" ")) > 1:
                                            intrefallele = map(int, refallele.split(" "))
                                            sortedrefallele = " ".join(str(allele) for allele in sorted(intrefallele))
                                        else:
                                            sortedrefallele = refallele
                                        # Populate self.mlstseqtype with the genome, best match to profile, # of matches
                                        # to the profile, gene, query allele(s), reference allele(s), and % identity
                                        if self.updateprofile:
                                            self.mlstseqtype[genome][sequencetype][sortedmatches][gene][
                                                str(self.bestdict[genome][gene]
                                                    .keys()[0])][sortedrefallele] = str(self.bestdict[genome][gene]
                                                                                        .values()[0])
                                        else:
                                            self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                                self.bestdict[genome][gene].keys()[0]] \
                                                = str(self.bestdict[genome][gene].values()[0])
                                            if sortedrefallele != self.bestdict[sample.name][gene].keys()[0]:
                                                mismatches.append(
                                                    ({gene: ('{} ({})'.format(self.bestdict[sample.name][gene]
                                                                              .keys()[0], sortedrefallele))}))
                                        if not self.updateprofile or self.analysistype == 'mlst':
                                            self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                                self.bestdict[genome][gene]
                                                    .keys()[0]] = str(self.bestdict[genome][gene].values()[0])
                                            sample[self.analysistype].mismatchestosequencetype = mismatches
                                            sample[self.analysistype].sequencetype = sequencetype
                                            sample[self.analysistype].matchestosequencetype = matches
                            # Add the new profile to the profile file (if the option is enabled)
                            if self.updateprofile and self.analysistype != 'mlst':
                                self.reprofiler(int(header), genome, sample)
                        elif sortedmatches == 0:
                            for gene in sample[self.analysistype].allelenames:
                                # Populate the results profile with negative values for sequence type and sorted matches
                                self.resultprofile[genome]['NA'][sortedmatches][gene]['NA'] = 0
                            # Add the new profile to the profile file (if the option is enabled)
                            if self.updateprofile:
                                self.reprofiler(int(header), genome, sample)
                            sample[self.analysistype].sequencetype = 'NA'
                            sample[self.analysistype].matchestosequencetype = 'NA'
                            sample[self.analysistype].mismatchestosequencetype = 'NA'
                        else:
                            sample[self.analysistype].matchestosequencetype = 'NA'
                            sample[self.analysistype].mismatchestosequencetype = 'NA'
                            sample[self.analysistype].sequencetype = 'NA'
                        dotter()
                else:
                    sample[self.analysistype].matchestosequencetype = 'NA'
                    sample[self.analysistype].mismatchestosequencetype = 'NA'
                    sample[self.analysistype].sequencetype = 'NA'

            else:
                sample[self.analysistype].matchestosequencetype = 'NA'
                sample[self.analysistype].mismatchestosequencetype = 'NA'
                sample[self.analysistype].sequencetype = 'NA'

    def reprofiler(self, header, genome, sample):
        """
        Creates and appends new profiles as required
        :param header:
        :param genome:
        :param sample:
        """
        # Iterate through mlstseqtype - it contains genomes with partial matches to current reference profiles
        # Reset :newprofile
        newprofile = ""
        # Find the last profile entry in the dictionary of profiles
        # Opens uses the command line tool 'tail' to look at the last line of the file (-1). This last line
        # is split on tabs, and only the first entry (the sequence type number) is captured
        if sample[self.analysistype].supplementalprofile != 'NA':
            if os.path.isfile(sample[self.analysistype].supplementalprofile):
                try:
                    lastentry = int(
                        subprocess.check_output(['tail', '-1', sample[self.analysistype].supplementalprofile])
                        .split("\t")[0]) + 1
                except ValueError:
                    lastentry = 1000000
            else:
                open(sample[self.analysistype].supplementalprofile, 'wb').close()
                lastentry = 1000000
            # As there can be multiple profiles in MLSTSeqType, this loop only needs to be performed once.
            seqcount = 0
            # Go through the sequence types
            try:
                sequencetype = self.mlstseqtype[genome].keys()[0]
            except IndexError:
                sequencetype = ''
                seqcount = 1
            # Only do this once
            if seqcount == 0:
                # Set the :newprofile string to start with the new profile name (e.g. 1000000_CFIA)
                newprofile = str(lastentry)
                # The number of matches to the reference profile
                nummatches = self.mlstseqtype[genome][sequencetype].keys()[0]
                # The genes in geneList - should be in the correct order
                for gene in sorted(sample[self.analysistype].allelenames):
                    # The allele for each gene in the query genome
                    allele = self.mlstseqtype[genome][sequencetype][nummatches][gene].keys()[0]
                    # Append the allele to newprofile
                    newprofile += '\t{}'.format(allele)
                    # Add the MLST results for the query genome as well as the new profile data
                    # to resultProfile
                    self.resultprofile[genome]['{}(new)'.format(str(lastentry))][header][gene][allele] = \
                        self.mlstseqtype[genome][sequencetype][nummatches][gene][allele].values()[0]
                seqcount += 1
                sample[self.analysistype].mismatchestosequencetype = 'NA'
                sample[self.analysistype].matchestosequencetype = header
            # Only perform the next loop if :newprofile exists
            if newprofile:
                # Open the profile file to append
                with open(sample[self.analysistype].supplementalprofile, 'ab') as appendfile:
                    # Append the new profile to the end of the profile file
                    appendfile.write('{}\n'.format(newprofile))
                # Re-run profiler with the updated files
                self.profiler()
        else:
            sample[self.analysistype].mismatchestosequencetype = 'NA'
            sample[self.analysistype].matchestosequencetype = 'NA'

    def reporter(self):
        """ Parse the results into a report"""
        # Initialise variables
        combinedrow = ''
        reportdirset = set()
        # Populate a set of all the report directories to use. A standard analysis will only have a single report
        # directory, while pipeline analyses will have as many report directories as there are assembled samples
        for sample in self.metadata:
            # Ignore samples that lack a populated reportdir attribute
            if sample[self.analysistype].reportdir != 'NA':
                make_path(sample[self.analysistype].reportdir)
                # Add to the set - I probably could have used a counter here, but I decided against it
                reportdirset.add(sample[self.analysistype].reportdir)
        # Create a report for each sample from :self.resultprofile
        for sample in self.metadata:
            if sample[self.analysistype].reportdir != 'NA':
                if type(sample[self.analysistype].allelenames) == list:
                    # Populate the header with the appropriate data, including all the genes in the list of targets
                    row = 'Strain,Genus,SequenceType,Matches,{},\n' \
                        .format(','.join(sorted(sample[self.analysistype].allelenames)))
                    # Set the sequence counter to 0. This will be used when a sample has multiple best sequence types.
                    # The name of the sample will not be written on subsequent rows in order to make the report clearer
                    seqcount = 0
                    # Iterate through the best sequence types for the sample (only occurs if update profile is disabled)
                    for seqtype in self.resultprofile[sample.name]:
                        """
                        {
                            "OLF15230-1_2015-SEQ-0783": {
                                "1000004_CFIA": {
                                    "7": {
                                        "dnaE": {
                                            "47": "100.00"
                                        },
                                        "dtdS": {
                                            "19": "100.00"
                                        },
                                        "gyrB": {
                                            "359": "100.00"
                                        },
                                        "pntA": {
                                            "50": "100.00"
                                        },
                                        "pyrC": {
                                            "143": "100.00"
                                        },
                                        "recA": {
                                            "31": "100.00"
                                        },
                                        "tnaA": {
                                            "26": "100.00"
                                        }
                                    }
                                }
                            }
                        }
                        """
                        # Becomes
                        """
                        Strain,SequenceType,Matches,dnaE,gyrB,recA,dtdS,pntA,pyrC,tnaA
                        OLF15230-1_2015-SEQ-0783,1000004_CFIA,7,26 (100.00%),359 (100.00%),31 (100.00%),50 (100.00%),
                            19 (100.00%),47 (100.00%),143 (100.00%)
                        """
                        sample[self.analysistype].sequencetype = seqtype
                        # The number of matches to the profile
                        matches = self.resultprofile[sample.name][seqtype].keys()[0]
                        # If this is the first of one or more sequence types, include the sample name
                        if seqcount == 0:
                            row += '{},{},{},{},'.format(sample.name, sample.general.referencegenus, seqtype, matches)
                        # Otherwise, skip the sample name
                        else:
                            row += ',,{},{},'.format(seqtype, matches)
                        # Iterate through all the genes present in the analyses for the sample
                        for gene in sorted(sample[self.analysistype].allelenames):
                            # refallele = self.profiledata[self.analysistype][seqtype][gene]
                            refallele = sample[self.analysistype].profiledata[seqtype][gene]
                            # Set the allele and percent id from the dictionary's keys and values, respectively
                            allele = self.resultprofile[sample.name][seqtype][matches][gene].keys()[0]
                            percentid = self.resultprofile[sample.name][seqtype][matches][gene].values()[0]
                            if refallele and refallele != allele:
                                if 0 < float(percentid) < 100:
                                    row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                else:
                                    row += '{} ({}),'.format(allele, refallele)
                            else:
                                # Add the allele and % id to the row (only add the percent identity if it is not 100%)
                                if 0 < float(percentid) < 100:
                                    row += '{} ({:.2f}%),'.format(allele, float(percentid))
                                else:
                                    row += '{},'.format(allele)
                            self.referenceprofile[sample.name][gene] = allele
                        # Add a newline
                        row += '\n'
                        # Increment the number of sequence types observed for the sample
                        seqcount += 1
                    combinedrow += row
                    # If the length of the # of report directories is greater than 1 (script is being run as part of
                    # the assembly pipeline) make a report for each sample
                    if self.pipeline:
                        # Open the report
                        with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                       self.analysistype), 'wb') as report:
                            # Write the row to the report
                            report.write(row)
                dotter()
            # Create the report folder
            make_path(self.reportpath)
            # Create the report containing all the data from all samples
            if self.pipeline:
                with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                        as combinedreport:
                    # Write the results to this report
                    combinedreport.write(combinedrow)
            else:
                with open('{}{}_{:}.csv'.format(self.reportpath, self.analysistype, time.strftime("%Y.%m.%d.%H.%M.%S")),
                          'wb') as combinedreport:
                    # Write the results to this report
                    combinedreport.write(combinedrow)
            # Remove the raw results csv
            [os.remove(rawresults) for rawresults in glob('{}*rawresults*'.format(self.reportpath))]

    def dumper(self):
        """Write :self.referenceprofile"""
        with open('{}{}_referenceprofile.json'.format(self.reportpath, self.analysistype, ), 'wb') as referenceprofile:
            referenceprofile.write(json.dumps(self.referenceprofile, sort_keys=True, indent=4, separators=(',', ': ')))

    def referencegenomefinder(self):
        """
        Finds the closest reference genome to the profile of interest
        """
        # Initialise dictionaries
        referencematch = defaultdict(make_dict)
        referencehits = defaultdict(make_dict)
        # Set the name of the reference profile file
        referencegenomeprofile = '{}rMLST_referenceprofile.json'.format(self.referenceprofilepath)
        # Open the reference profile and load the profile into memory
        with open(referencegenomeprofile) as referencefile:
            referencetypes = json.load(referencefile)
        # Iterate through the samples
        for sample in self.metadata:
            if sample[self.analysistype].reportdir != 'NA':
                # Iterate through the reference genomes in the profile
                for genome in referencetypes:
                    # Initialise the number of identical alleles between the assembly of interest and the
                    # reference genome to 0
                    referencehits[sample.name][genome] = 0
                    # Iterate through all the genes in the analysis
                    for gene in self.bestdict[sample.name]:
                        # If the alleles match between the assembly of interest and the reference genome, increment
                        # the number of matches by 1
                        if self.bestdict[sample.name][gene].keys()[0] == referencetypes[genome][gene]:
                            referencematch[sample.name][genome][gene] = 1
                            referencehits[sample.name][genome] += 1
                        else:
                            referencematch[sample.name][genome][gene] = 0
        #
        for sample in self.metadata:
            if sample[self.analysistype].reportdir != 'NA':
                # Get the best number of matches
                # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
                try:
                    sortedmatches = sorted(referencehits[sample.name].items(),
                                           key=operator.itemgetter(1), reverse=True)[0]
                except IndexError:
                    sortedmatches = (0, 0)
                # If there are fewer matches than the total number of genes in the typing scheme
                if 0 < int(sortedmatches[1]) < len(sample[self.analysistype].allelenames):
                    mismatches = []
                    # Iterate through the gene in the analysis
                    for gene, allele in referencetypes[sortedmatches[0]].iteritems():
                        # Populate :self.referencegenome with the genome name, best reference match, number of matches,
                        # gene, query allele(s), and percent identity
                        percentidentity = '{:.2f}'.format(self.bestdict[sample.name][gene].values()[0])
                        self.referencegenome[sample.name][sortedmatches[0]][sortedmatches[1]][gene][self.bestdict[
                            sample.name][gene].keys()[0]] = percentidentity
                        if self.bestdict[sample.name][gene].keys()[0] != allele:
                            sample[self.analysistype].referencegenome = sortedmatches[0]
                            sample.general.referencegenus = sortedmatches[0].split('_')[0]
                            sample[self.analysistype].referencegenomepath = '{}{}.fa' \
                                .format(self.referenceprofilepath, sortedmatches[0])
                            sample[self.analysistype].matchestoreferencegenome = sortedmatches[1]
                            mismatches.append(({gene: ('{} ({})'.format(self.bestdict[sample.name][gene]
                                                                        .keys()[0], allele))}))
                        sample[self.analysistype].mismatchestoreferencegenome = sorted(mismatches)
                elif sortedmatches == 0:
                    for gene in sample[self.analysistype].allelenames:
                        # Populate the profile of results with 'negative' values for sequence type and sorted matches
                        self.referencegenome[sample.name][sortedmatches[0]][0][gene]['NA'] = 0
                        sample[self.analysistype].referencegenome = 'NA'
                        sample.general.referencegenus = 'NA'
                        sample[self.analysistype].referencegenomepath = 'NA'
                        sample[self.analysistype].matchestoreferencegenome = 0
                        sample[self.analysistype].mismatchestoreferencegenome = [0]
                # Otherwise, the query profile matches the reference profile
                else:
                    for gene in referencetypes[sortedmatches[0]]:
                        # Populate self.referencegenome as above
                        self.referencegenome[sample.name][sortedmatches[0]][sortedmatches[1]][gene][self.bestdict[
                            sample.name][gene].keys()[0]] = '{:.2f}'.format(self.bestdict[
                                                                                sample.name][gene].values()[0])
                        sample[self.analysistype].referencegenome = sortedmatches[0]
                        sample[self.analysistype].referencegenomepath = '{}{}.fa' \
                            .format(self.referenceprofilepath, sortedmatches[0])
                        sample.general.referencegenus = sortedmatches[0].split('_')[0]
                        sample[self.analysistype].matchestoreferencegenome = sortedmatches[1]
                        sample[self.analysistype].mismatchestoreferencegenome = [0]
        # Print the results to file
        make_path(self.reportpath)
        with open('{}referencegenomes.csv'.format(self.reportpath), 'wb') as referencegenomereport:
            row = 'Strain,referencegenome,numberofmatches\n'
            for sample in self.metadata:
                if sample[self.analysistype].reportdir != 'NA':
                    row += '{},{},{}\n'.format(sample.name, sample[self.analysistype].referencegenome,
                                               sample[self.analysistype].matchestoreferencegenome)
            referencegenomereport.write(row)
            dotter()

    def __init__(self, inputobject):
        from threading import Lock
        import multiprocessing
        self.path = inputobject.path
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = inputobject.start
        self.analysistype = inputobject.analysistype
        self.allelefolders = set()
        self.fnull = open(os.devnull, 'w')  # define /dev/null
        self.lock = Lock()
        self.updateallele = inputobject.updateallele
        self.updateprofile = inputobject.updateprofile
        self.updatedb = []
        self.reportpath = inputobject.reportdir
        self.datadump = inputobject.datadump
        self.bestreferencegenome = inputobject.bestreferencegenome
        self.pipeline = inputobject.pipeline
        self.referencefilepath = inputobject.referencefilepath
        self.referenceprofilepath = inputobject.referenceprofilepath
        # Fields used for custom outfmt 6 BLAST output:
        # "6 qseqid sseqid positive mismatch gaps evalue bitscore slen length"
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length']
        self.cpus = int(multiprocessing.cpu_count())
        # Declare queues, and dictionaries
        self.dqueue = Queue(maxsize=self.cpus)
        self.blastqueue = Queue(maxsize=self.cpus)
        self.blastdict = {}
        self.blastresults = defaultdict(make_dict)
        self.plusdict = defaultdict(make_dict)
        self.bestdict = defaultdict(make_dict)
        self.bestmatch = defaultdict(int)
        self.mlstseqtype = defaultdict(make_dict)
        self.resultprofile = defaultdict(make_dict)
        # self.profiledata = defaultdict(make_dict)
        self.referenceprofile = defaultdict(make_dict)
        self.referencegenome = defaultdict(make_dict)
        # Run the MLST analyses
        self.mlst()


def allelesplitter(allelenames):
    # Multiple try-excepts. Maybe overly complicated, but I couldn't get it work any other way
    # This (hopefully) accounts for all the possible naming schemes for the alleles
    try:  # no split - just allele numbers e.g. >12
        match = re.search(r"(>\d+)", allelenames)
        allelenumber = str(match.group()).split(">")[1]
        alleleprenumber = ""
    except (IndexError, AttributeError):
        try:  # split on "_" e.g. >AROC_12
            # allelenumber is the number of the allele(!). It should be different for each allele
            allelenumber = allelenames.split("_")[1]
            # alleleprenumber is anything before the allele number. It should be the same for each allele
            alleleprenumber = allelenames.split("_")[0]
        except IndexError:
            try:  # split on "-" e.g. >AROC-12
                allelenumber = allelenames.split("-")[1]
                alleleprenumber = allelenames.split("-")[0]
            except IndexError:
                try:  # split on " " e.g. >AROC 12
                    allelenumber = allelenames.split(" ")[1]
                    alleleprenumber = allelenames.split(" ")[0]
                except IndexError:
                    try:  # split on change from letters to numbers e.g. >AROC12
                        match = re.search(r"(>[A-z/a-z]+)(\d+)", allelenames)
                        allelenumber = match.groups()[1]
                        alleleprenumber = match.groups()[0]
                    except (IndexError, AttributeError):
                        allelenumber = allelenames
                        alleleprenumber = allelenames
    # Return the variables
    return int(allelenumber), alleleprenumber


def blastdatabaseclearer(genepath):
    """
    Due to the nature of the program updating allele files, it's not desirable to use previously generated databases.
    Additionally, with the use of these files by multiple programs, there is an issue. This script makes database files
    as follows: aroC.fasta becomes aroC.nhr, etc. The current SPAdes assembly pipeline would take that same .fasta file
    and create aroC.fasta.nhr. Deleting database files prevents issues with glob including database files.
    :param genepath: path to folder containing the MLST target genes
    """
    # Get all the .nhr, .nin, .nsq files
    databaselist = glob('{}/*.n*'.format(genepath))
    # And delete them
    for allele in databaselist:
        os.remove(allele)


def combinealleles(start, allelepath, alleles):
    printtime('Creating combined rMLST allele file', start)
    with open('{}/rMLST_combined.fasta'.format(allelepath), 'wb') as combinedfile:
        # Open each allele file
        for allele in sorted(alleles):
            # with open(allele, 'rU') as fasta:
            for record in SeqIO.parse(open(allele, "rU"), "fasta"):
                # Extract the sequence record from each entry in the multifasta
                # Replace and dashes in the record.id with underscores
                record.id = record.id.replace('-', '_')
                # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                # with gaps
                # noinspection PyProtectedMember
                record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                # Clear the name and description attributes of the record
                record.name = ''
                record.description = ''
                # Write each record to the combined file
                SeqIO.write(record, combinedfile, 'fasta')


def schemedate(lastfolder):
    from datetime import date
    try:
        # Extract the folder name (date) from the path/name
        lastupdate = os.path.split(lastfolder)[1]
    except AttributeError:
        lastupdate = '2000-01-01'
    try:
        # Calculate the size of the folder by adding the sizes of all the files within the folder together
        foldersize = sum(os.path.getsize('{}/{}'.format(lastfolder, f)) for f in os.listdir(lastfolder)
                         if os.path.isfile('{}/{}'.format(lastfolder, f)))
    except TypeError:
        foldersize = 0
    # Try to figure out the year, month, and day from the folder name
    try:
        (year, month, day) = lastupdate.split("-")
        # Create a date object variable with the year, month, and day
        d0 = date(int(year), int(month), int(day))
    except ValueError:
        # Set an arbitrary date in the past to force an update
        d0 = date(2000, 01, 01)
    # Create a date object with the current date
    d1 = date(int(time.strftime("%Y")), int(time.strftime("%m")), int(time.strftime("%d")))
    # Subtract the last update date from the current date
    delta = d1 - d0

    return delta, foldersize, d1


def getrmlsthelper(referencefilepath, update, start):
    """
    Makes a system call to rest_auth.pl, a Perl script modified from
    https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
    And downloads the most up-to-date rMLST profile and alleles
    """
    from subprocess import call
    analysistype = 'rMLST'
    # Folders are named based on the download date e.g 2016-04-26
    # Find all folders (with the trailing / in the glob search) and remove the trailing /
    lastfolder = sorted(glob('{}{}/*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
    delta, foldersize, d1 = schemedate(lastfolder)
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Set the path/name of the folder to contain the new alleles and profile
    newfolder = '{}{}/{}'.format(referencefilepath, analysistype, d1)
    # System call
    rmlstupdatecall = 'cd {} && perl {}/rest_auth.pl -a {}/secret.txt'.format(newfolder, homepath, homepath)
    if update:
        if delta.days > 7 or foldersize < 100:
            printtime("Last update of rMLST profile and alleles was {} days ago. Updating".format(str(delta.days)),
                      start)
            # Create the path
            make_path(newfolder)
            # Copy over the access token to be used in the authentication
            shutil.copyfile('{}/access_token'.format(homepath), '{}/access_token'.format(newfolder))
            # Run rest_auth.pl
            call(rmlstupdatecall, shell=True)
            # Get the new alleles into a list, and create the combinedAlleles file
            alleles = glob('{}/*.tfa'.format(newfolder))
            combinealleles(start, newfolder, alleles)
        # If the profile and alleles are up-to-date, set :newfolder to :lastfolder
        else:
            newfolder = lastfolder
        # Ensure that the profile/alleles updated successfully
        # Calculate the size of the folder by adding the sizes of all the files within the folder together
        newfoldersize = sum(os.path.getsize('{}/{}'.format(newfolder, f)) for f in os.listdir(newfolder)
                            if os.path.isfile('{}/{}'.format(newfolder, f)))
        # If the profile/allele failed, remove the folder, and use the most recent update
        if newfoldersize < 100:
            shutil.rmtree(newfolder)
            newfolder = sorted(glob('{}{}/*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
    # Don't update the profile/alleles if not requested
    else:
        newfolder = lastfolder
    # Return the system call and the folder containing the profile and alleles
    return rmlstupdatecall, newfolder


def getmlsthelper(referencefilepath, start, organism, update):
    """Prepares to run the getmlst.py script provided in SRST2"""
    from accessoryFunctions import GenObject
    # Initialise a set to for the organism(s) for which new alleles and profiles are desired
    organismset = set()
    # Allow for Shigella to use the Escherichia MLST profile/alleles
    organism = organism if organism != 'Shigella' else 'Escherichia'
    # As there are multiple profiles for certain organisms, this dictionary has the schemes I use as values
    organismdictionary = {'Escherichia': 'Escherichia coli#1',
                          'Shigella': 'Escherichia coli#1',
                          'Vibrio': 'Vibrio parahaemolyticus',
                          'Campylobacter': 'Campylobacter jejuni',
                          'Listeria': 'Listeria monocytogenes',
                          'Bacillus': 'Bacillus cereus'}
    # Allow for a genus not in the dictionary being specified
    try:
        organismset.add(organismdictionary[organism])
    except KeyError:
        # Add the organism to the set
        organismset.add(organism)
    for scheme in organismset:
        organismpath = os.path.join(referencefilepath, 'MLST', organism)
        # Find all folders (with the trailing / in the glob search) and remove the trailing /
        try:
            lastfolder = sorted(glob('{}/*/'.format(organismpath)))[-1].rstrip('/')
        except IndexError:
            lastfolder = []
        # Run the method to determine the most recent folder, and how recently it was updated
        delta, foldersize, d1 = schemedate(lastfolder)
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = '{}/{}'.format(organismpath, d1)
        if update:
            if delta.days > 7 or foldersize < 100:
                printtime('Downloading {} MLST scheme from pubmlst.org'.format(organism), start)
                # Create the object to store the argument attributes to feed to getmlst
                getmlstargs = GenObject()
                getmlstargs.species = scheme
                getmlstargs.repository_url = 'http://pubmlst.org/data/dbases.xml'
                getmlstargs.force_scheme_name = False
                getmlstargs.path = newfolder
                # Create the path to store the downloaded
                make_path(getmlstargs.path)
                getmlst.main(getmlstargs)
                # Even if there is an issue contacting the database, files are created, however, they are populated
                # with XML strings indicating that the download failed
                # Read the first character in the file
                try:
                    profilestart = open(glob('{}/*.txt'.format(newfolder))[0]).readline()
                except IndexError:
                    profilestart = []
                # If it is a <, then the download failed
                if not profilestart or profilestart[0] == '<':
                    # Delete the folder, and use the previous definitions instead
                    shutil.rmtree(newfolder)
                    newfolder = lastfolder
            # If the profile and alleles are up-to-date, set :newfolder to :lastfolder
            else:
                newfolder = lastfolder
        # If update isn't specified, don't update
        else:
            newfolder = lastfolder
            # Ensure that the profile/alleles updated successfully
            # Calculate the size of the folder by adding the sizes of all the files within the folder together
        try:
            newfoldersize = sum(os.path.getsize('{}/{}'.format(newfolder, f)) for f in os.listdir(newfolder)
                                if os.path.isfile('{}/{}'.format(newfolder, f)))
        except (OSError, TypeError):
            newfoldersize = 100
        # If the profile/allele failed, remove the folder, and use the most recent update
        if newfoldersize < 100:
            shutil.rmtree(newfolder)
            try:
                newfolder = sorted(glob('{}/*/'.format(organismpath)))[-1].rstrip('/')
            except IndexError:
                newfolder = organismpath
        # Return the name/path of the allele-containing folder
        return newfolder

if __name__ == '__main__':
    class Parser(object):

        def strainer(self):
            from accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = sorted(glob('{}*.fa*'.format(self.sequencepath))) if self.sequencepath \
                else sorted(glob('{}sequences/*.fa*'.format(self.path)))
            # Populate the metadata object. This object will be populated to mirror the objects created in the
            # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
            # of a pipeline
            assert self.strains, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your sequence path is correct'.format(self.sequencepath)
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
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def organismchooser(self):
            """Find the alleles and profiles based on inputs"""
            if self.analysistype == 'MLST':
                if not self.allelepath:
                    # If the name of the organism to analyse was provided
                    assert self.organism, 'Need to provide either a path to the alleles or an organism name'
                    # If the -g flag was included, download the appropriate MLST scheme for the organism
                    # if self.getmlst and self.organism:
                    # referencefilepath, start, scheme, path, organism
                    # getmlsthelper(self.r)
                    self.allelepath = '{}{}'.format(self.path, self.organism)
                    assert os.path.isdir(self.allelepath), 'Cannot find {}. Please ensure that the folder exists, or ' \
                                                           'use the -g option to download the {} MLST scheme' \
                        .format(self.allelepath, self.organism)
                # Tries to get the organism name for the folder containing the alleles
                self.organism = self.organism if self.organism else os.path.split(self.allelepath)[-1]
                if self.cleardatabases:
                    # Remove and previously created blast database files (if desired)
                    blastdatabaseclearer(self.allelepath)
                # Create lists of the alleles, combined alleles, and the profile
                self.alleles = glob('{}/*.tfa'.format(self.allelepath)) if glob('{}/*.tfa'.format(self.allelepath)) \
                    else glob('{}/*.fas'.format(self.allelepath))
                self.combinedalleles = glob('{}/*.fasta'.format(self.allelepath))
                # Get the .txt profile file name and path into a variable
                self.profile = glob('{}/*.txt'.format(self.allelepath))
            # rMLST analyses are slightly different; alleles cannot be downloaded in the same fashion as MLST alleles,
            # and the alleles have different file extensions
            else:
                if not self.allelepath:
                    self.allelepath = '{}rMLST'.format(self.path)
                    # If the name of the organism to analyse was provided
                    assert os.path.isdir(self.allelepath), 'Cannot find directory containing rMLST alleles and ' \
                                                           'profile in {}'.format(self.path)
                # Create lists of the alleles, combined alleles, and the profile
                self.alleles = glob('{}*.fas'.format(self.allelepath))
                self.combinedalleles = glob('{}*.fasta'.format(self.allelepath))
                # If the combined alleles files doesn't exist
                size = 0
                if self.combinedalleles:
                    size = os.stat(self.combinedalleles[0]).st_size
                if not self.combinedalleles or size == 0:
                    # Open the combined allele file to write
                    combinealleles(self.start, self.allelepath, self.alleles)
                # Set the combined alleles file name and path
                self.combinedalleles = glob('{}/*.fasta'.format(self.allelepath))
                # Get the .txt profile file name and path into a variable
                self.profile = glob('{}/*.txt'.format(self.allelepath))
            for sample in self.samples:
                sample[self.analysistype].alleles = [os.path.split(x)[1].split('.')[0] for x in self.alleles]
                sample[self.analysistype].allelenames = [os.path.split(x)[1].split('.')[0] for x in self.alleles]
                sample[self.analysistype].alleledir = self.allelepath
                sample[self.analysistype].profile = self.profile
                sample[self.analysistype].analysistype = self.scheme
                sample[self.analysistype].reportdir = self.reportpath
                sample[self.analysistype].organism = self.organism
                sample[self.analysistype].combinedalleles = self.combinedalleles

        def __init__(self):
            from argparse import ArgumentParser
            parser = ArgumentParser(description='Performs blast analyses to determine presence of alleles in a genome '
                                                'query, and types genome based on typing profile. Adds novel alleles '
                                                'and profiles to the appropriate files. '
                                                'Example command: '
                                                'python mMLST.py'
                                                '-p /home/git/MLST  '
                                                '-s /home/git/MLST/sequences '
                                                '-O /home/git/MLST/Organism '
                                                '-o Vibrio '
                                                '-S MLST')
            parser.add_argument('-p', '--path',
                                default=os.getcwd(),
                                help='Specify path for custom folder locations. If you don\'t supply additional paths'
                                     'e.g. sequencepath, allelepath, or organismpath, then the program will look for '
                                     'MLST files in .../path/Organism, and the query sequences in ../path/sequences. '
                                     'If you don\'t input a path, then the current working directory will be used.')
            parser.add_argument('-c', '--cutoff',
                                default=98,
                                help='The percent identity cutoff value for BLAST matches. Default is 98%)')
            parser.add_argument('-s', '--sequencepath',
                                help='The location of the query sequence files')
            parser.add_argument('-a', '--alleleprofilepath',
                                help='The path of the folder containing the allele files, and the profile file '
                                     'e.g. /folder/path/Organism/Vibrio/cgMLST'
                                     'Please note the requirements for the profile database in the readme')
            parser.add_argument('-O', '--organismpath',
                                help='The path of the folder containing the organism folders e.g. folder/path/Organism')
            parser.add_argument('-o', '--organism',
                                help='The name of the organism you wish to type. Must match the folder name containing '
                                     'the schemes e.g. Salmonella or "Clostridium botulinum" (note the quotes')
            parser.add_argument('-S', '--scheme',
                                help='The scheme you wish to use. Must match the folder name containing the scheme e.g.'
                                     ' cgMLST. Furthermore, this folder must contain two folders: "alleles" and '
                                     '"profile". The alleles folder contains the allele files in .fasta format, and the'
                                     ' profile folder contains the profile in .txt format. Please note the requirements'
                                     ' for the profile in the readme')
            parser.add_argument('-u', '--updateprofilefalse',
                                action='store_true',
                                help='By default, the program does not create new sequence profiles and appends '
                                     'these profiles to the profile file. Including this flag enables this '
                                     'functionality.')
            parser.add_argument('-U', '--updateallelefalse',
                                action='store_true',
                                help='By default, the program does not create new alleles and appends these '
                                     'alleles to the appropriate file. Including this flag enables this functionality')
            parser.add_argument('-r', '--reportdirectory',
                                default='{}/reports'.format(os.getcwd()),
                                help='Path to store the reports defaults to os.getcwd()/reports')
            parser.add_argument('-d', '--dumpdata',
                                action='store_true',
                                help='Optionally dump :self.resultprofile dictionary to file. Useful when creating a '
                                     'reference database against which novel sequences can be compared. '
                                     'The .json file will be placed in the reports folder ')
            parser.add_argument('-g', '--getmlst',
                                action='store_true',
                                help='Optionally get the newest profile and alleles for your analysis from pubmlst.org')
            parser.add_argument('-b', '--bestreferencegenome',
                                action='store_true',
                                help='Optionally find the refseq genome with the largest number of rMLST alleles in '
                                     'common with the strain of interest')
            parser.add_argument('-R', '--referenceprofile',
                                help='Path with reference genomes and a referenceprofile.json file (created with the'
                                     '-d flag. The genomes in the .json profile must be the genomes in the folder')
            parser.add_argument('-t', '--type',
                                default='MLST',
                                help='Specify the analysis type (MLST or rMLST now. More types of analysis may become'
                                     'available in the future')
            parser.add_argument('-C', '--clearblastdatabases',
                                action='store_true',
                                help='By default, the BLAST database for your analysis are not deleted prior the '
                                     'analyses. Potentially, the most up-to-date allele definitions will not be used. '
                                     'Use the -C flag to enable the deletion of the databases')

            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.reportpath = os.path.join(args.reportdirectory, '')
            self.cutoff = float(args.cutoff)
            self.sequencepath = os.path.join(args.sequencepath, '') if args.sequencepath else '{}sequences/' \
                .format(self.path) if os.path.isdir('{}sequences'.format(self.path)) else self.path
            self.allelepath = os.path.join(args.alleleprofilepath, '') if args.alleleprofilepath else ''
            self.organismpath = os.path.join(args.organismpath, '') if args.organismpath else ''
            self.referenceprofilepath = os.path.join(args.referenceprofile, '') if args.referenceprofile else \
                '{}referenceGenomes/'.format(self.path)
            if args.bestreferencegenome:
                assert os.path.isdir(self.referenceprofilepath), 'Cannot find {}. Please double check that you ' \
                                                                 'provided the proper path to the reference profile ' \
                                                                 'folder'.format(self.referenceprofilepath)
            self.scheme = args.scheme if args.scheme else args.type
            self.organism = args.organism
            self.updateprofile = args.updateprofilefalse
            self.updateallele = args.updateallelefalse
            self.datadump = args.dumpdata
            self.getmlst = args.getmlst
            self.bestreferencegenome = args.bestreferencegenome
            self.analysistype = args.type
            self.cleardatabases = args.clearblastdatabases

            # Initialise variables
            self.genepath = ''
            self.alleles = ''
            self.combinedalleles = ''
            self.profile = ''
            self.strains = []
            self.samples = []
            self.start = time.time()
            # Get a list of the sequence files
            self.strainer()
            self.organismchooser()


    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.start = self.runmetadata.start
            self.analysistype = self.runmetadata.scheme if self.runmetadata.scheme else self.runmetadata.analysistype
            self.alleles = self.runmetadata.alleles
            self.profile = self.runmetadata.profile
            # self.supplementalprofile = False
            self.cutoff = self.runmetadata.cutoff
            self.updateallele = self.runmetadata.updateallele
            self.updateprofile = self.runmetadata.updateprofile
            self.reportdir = self.runmetadata.reportpath
            self.datadump = self.runmetadata.datadump
            self.getmlst = self.runmetadata.getmlst
            self.bestreferencegenome = self.runmetadata.bestreferencegenome
            self.referenceprofilepath = self.runmetadata.referenceprofilepath
            self.referencefilepath = self.runmetadata.allelepath
            self.pipeline = False
            # Run the analyses
            MLST(self)

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        """
        Determine whether it is required to run the MLST analyses
        """
        # Initialise a variable to store whether the analyses need to be performed
        analyse = list()
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                try:
                    # Try to open the final report from the analyses. If it exists, then the analyses don't need to be
                    # performed again.
                    if os.path.isfile('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                           self.analysistype)):
                        # Run the allele updater method
                        updatecall, allelefolder = getrmlsthelper(self.referencefilepath, self.updatedatabases,
                                                                  self.start)
                        # Alleles have a .tfa extension
                        self.alleles = glob('{}/*.tfa'.format(allelefolder))
                        sample[self.analysistype].alleles = self.alleles
                        sample[self.analysistype].allelenames = [os.path.split(x)[1].split('.')[0] for x in
                                                                 self.alleles]
                        # The analyses have already been successfully completed
                        analyse.append(False)
                    # Otherwise run the analyses
                    else:
                        self.populator(sample)
                        analyse.append(True)
                # If the attribute doesn't exist, then the analyses haven't been performed yet.
                except (KeyError, AttributeError):
                    self.populator(sample)
                    analyse.append(True)
            else:
                self.populator(sample)
                analyse.append(True)
        # Only run the analyses if they have not completed successfully before
        if any(analyse):
            # Run the MLST analyses
            MLST(self)

    def populator(self, sample):
        """
        Populates objects with the necessary attributes
        :param sample: sample object
        """
        from accessoryFunctions import GenObject
        if sample.general.bestassemblyfile != 'NA':
            profile = ''
            setattr(sample, self.analysistype, GenObject())
            sample[self.analysistype].analyse = True
            if self.analysistype.lower() == 'rmlst':
                # Run the allele updater method
                updatecall, allelefolder = getrmlsthelper(self.referencefilepath, self.updatedatabases, self.start)
                # Alleles have a .tfa extension
                self.alleles = glob('{}/*.tfa'.format(allelefolder))
                # Get the profile file into a list
                profile = glob('{}/*.txt'.format(allelefolder))
                self.supplementalprofile = '{}rMLST/OLC_rMLST_profiles.txt'.format(self.referencefilepath)
                self.combinedalleles = glob('{}/*.fasta'.format(allelefolder))
                # Set the metadata file appropriately
                sample[self.analysistype].alleledir = allelefolder
                sample[self.analysistype].updatecall = updatecall
            else:
                # Use the getmlsthelper module to download databases as required
                schemefolder = getmlsthelper(self.referencefilepath, self.start, sample.general.referencegenus,
                                             self.updatedatabases)
                # If there is no database folder, do not perform the MLST analyses on this sample
                if not schemefolder:
                    sample[self.analysistype].analyse = False
                    # Set the metadata file appropriately
                    sample[self.analysistype].alleles = 'NA'
                    sample[self.analysistype].allelenames = 'NA'
                    sample[self.analysistype].profile = 'NA'
                    sample[self.analysistype].analysistype = 'NA'
                    sample[self.analysistype].reportdir = 'NA'
                    sample[self.analysistype].combinedalleles = 'NA'
                    sample[self.analysistype].supplementalprofile = 'NA'
                    sample[self.analysistype].alleledir = 'NA'
                else:
                    self.alleles = glob('{}/*.tfa'.format(schemefolder))
                    profile = glob('{}/*.txt'.format(schemefolder))
                    self.combinedalleles = glob('{}/*.fasta'.format(schemefolder))
                    sample[self.analysistype].alleledir = schemefolder
            # Only perform analyses on samples with target databases
            if sample[self.analysistype].analyse:
                sample[self.analysistype].alleles = self.alleles
                sample[self.analysistype].allelenames = [os.path.split(x)[1].split('.')[0] for x in self.alleles]
                sample[self.analysistype].profile = profile if profile else 'NA'
                sample[self.analysistype].analysistype = self.analysistype
                sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                      self.analysistype)
                sample[self.analysistype].combinedalleles = self.combinedalleles
                sample[self.analysistype].supplementalprofile = self.supplementalprofile \
                    if self.supplementalprofile else 'NA'

        else:
            setattr(sample, self.analysistype, GenObject())
            # Set the metadata file appropriately
            sample[self.analysistype].alleledir = 'NA'
            sample[self.analysistype].alleles = 'NA'
            sample[self.analysistype].allelenames = 'NA'
            sample[self.analysistype].profile = 'NA'
            sample[self.analysistype].analysistype = 'NA'
            sample[self.analysistype].reportdir = 'NA'
            sample[self.analysistype].combinedalleles = 'NA'
            sample[self.analysistype].supplementalprofile = 'NA'

    def __init__(self, inputobject, analysistype):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.referencefilepath = inputobject.reffilepath
        self.reportdir = '{}/'.format(inputobject.reportpath)
        self.alleles = ''
        self.profile = ''
        self.supplementalprofile = ''
        self.combinedalleles = ''
        self.cutoff = 99
        self.updatedatabases = inputobject.updatedatabases
        self.updateallele = True
        self.updateprofile = True
        self.datadump = False
        self.bestreferencegenome = True
        self.pipeline = True
        self.updatermlst = False
        self.referenceprofilepath = '{}referenceGenomes/'.format(self.referencefilepath)
        # Get the alleles and profile into the metadata
        self.strainer()
