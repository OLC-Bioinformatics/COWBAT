#! /usr/env/python
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
from cStringIO import StringIO
from glob import glob
import subprocess, os, time, shlex, re, threading, json, operator

from accessoryFunctions import dotter, make_path, make_dict, printtime

__author__ = 'akoziol, mikeknowles'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
Revised with speed improvements
"""


# TODO reportPath argument
# TODO account for 0 matches


class MLST(object):
    def mlst(self):
        # Get the MLST profiles into a dictionary for each sample
        self.profiler()
        # Make blast databases (if necessary)
        self.makedbthreads(self.allelefolders)
        # Run the blast analyses
        self.blastnthreads()
        # Determine sequence types from the analyses
        self.sequencetyper()
        # print json.dumps(self.resultprofile, sort_keys=True, indent=4, separators=(',', ': '))
        # Create reports
        self.reporter()

    def profiler(self):
        """Creates a dictionary from the profile scheme(s)"""
        # Initialise the dictionary
        profiledata = defaultdict(make_dict)
        for sample in self.metadata:
            genelist = []
            if sample.mlst.analysistype.keys()[0] == self.analysistype:
                if sample.mlst.profile[self.analysistype] != 'NA':
                    # The gene names are present in the first line of the profile file
                    # Get the MLST profiles for each sequence type
                    with open(sample.mlst.profile[self.analysistype][0]) as profile:
                        # Files have to be in tab-delimited format
                        header = profile.readline().rstrip().split("\t")
                        # As certain typing schemes are not in any discernible order, using a naturally ordered list
                        # instead of a dictionary to store the names is a good idea
                        for gene in header:
                            # The first column must have "ST" in the header
                            if 'ST' not in gene:
                                genelist.append(gene)
                                dotter()
                        for line in profile:
                            # Grab the name of the last profile
                            # mlstcount will used to associate the gene name in header to the allele (e.g. adk 12)
                            mlstcount = 1
                            # Don't need to get the header information again
                            if 'ST' not in line:
                                # len(header) will be the same length as the data in line
                                while mlstcount < len(header):
                                    # Remove newlines and split on tabs
                                    data = line.rstrip().split('\t')
                                    # Populate profileData with the sequence type, gene name, and the allele number
                                    profiledata[data[0]][header[mlstcount]] = data[mlstcount]
                                    # Increment the count
                                    mlstcount += 1
                    # Populate the metadata with the profile data
                    sample.mlst.profiledata = {self.analysistype: profiledata}
                    # Add the allele directory to a list of directories used in this analysis
                    self.allelefolders.add(sample.mlst.alleledir[self.analysistype])
                    # Add the list of genes in the analysis to each sample
                    sample.mlst.genelist = {self.analysistype: genelist}

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
            allelefiles = glob('{}*fa*'.format(alleledir))
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
            fnull = open(os.devnull, 'w')  # define /dev/null
            if not os.path.isfile(str(nhr)):  # if check for already existing dbs
                # Create the databases
                subprocess.call(shlex.split('makeblastdb -in {} -dbtype nucl -out {}'.format(fastapath, db)),
                                stdout=fnull, stderr=fnull)
            dotter()
            self.dqueue.task_done()  # signals to dqueue job is done

    def blastnthreads(self):
        """Setup and create  threads for blastn and xml path"""
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                for i in range(len(sample.mlst.alleles[self.analysistype])):
                    threads = Thread(target=self.runblast, args=())
                    threads.setDaemon(True)
                    threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                for allele in sample.mlst.alleles[self.analysistype]:
                    # Add each fasta/allele file combination to the threads
                    self.blastqueue.put((sample.general.bestassemblyfile, allele))
        # Join the threads
        self.blastqueue.join()

    def runblast(self):
        while True:  # while daemon
            (assembly, allele) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            genename = os.path.split(allele)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db,
            # a mild evalue of 0.1, and XML formatted output
            # Removed perc_identity=percentIdentity from the call, as this allows more flexibility for parsing files
            # with different cutoff values - if I want to re-analyse a search with a lower cutoff, I won't have to
            # re-perform the BLAST search each time
            db = allele.split(".")[0]
            blastn = NcbiblastnCommandline(query=assembly, db=db, evalue=0.1, outfmt=5)
            # Note that there is no output file specified -  the search results are currently stored in stdout
            stdout, stderr = blastn()
            # Search stdout for matches - if the term Hsp appears (the .find function will NOT
            # return -1), a match has been found, and stdout is written to file
            if stdout.find('Hsp') != -1:
                blast_handle = StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
                self.blastparse(blast_handle, genome, genename, allele)  # parse the data already in memory
            # If there are no records, populate the dictionary with "N" and a 0% identity
            else:
                # _0
                self.plusdict[genome]['{}'.format(genename)]["N"] = 0
            # dotter()
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparse(self, blast_handle, genome, gene, genepath):
        """
        Parses BLAST results, and populates a dictionary with the results
        :param blast_handle:
        :param genome:
        :param gene
        :param genepath
        """
        # global plusdict
        # global profileData
        # global updateallele
        snpdict = {}
        datadict = {}
        records = NCBIXML.parse(blast_handle)  # Open record from memory-mapped file
        dotter()
        # Split the extension from the genome name
        numhsp = sum(line.count('<Hsp>') for line in iter(blast_handle.readline, ""))
        if numhsp >= 1:
            # Since we scanned through result_handle looking for HSPs, the position of the read/write pointer
            # within the file is at the end. To reset it back to the beginning, .seek(0) is used
            blast_handle.seek(0)
            # Iterate through each record in the memory-mapped xml file
            for record in records:  # This process is just to retrieve HSPs from xml files
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        # Calculate the percent identity
                        percentidentity = '{:.2f}'.format(float(float(hsp.identities) / float(alignment.length) * 100))
                        allele = str(alignment.title.split(" ")[-1])
                        # If the results are 100% identical to the reference allele, add them to the dictionary
                        if hsp.identities >= alignment.length:
                            # Clears out any "N" values in the dictionary
                            if "N" in self.plusdict[genome][gene]:
                                self.plusdict[genome][gene].clear()
                            self.plusdict[genome][gene][allele] = percentidentity
                            # As the blast results files are not sorted by percent identity, and, at least for rMLST
                            # genes, not all genes are the same length, a match can occur after lots of close matches
                            snpdict.clear()
                        # Optionally process results if the percent identity is cutoff-99.9%
                        # and there is no 100% match in the dictionary
                        elif gene not in self.plusdict[genome] and not snpdict \
                                and hsp.identities >= alignment.length * self.cutoff:
                            if self.updateallele:
                                # Puts the HSP in the correct order -  hits to the negative strand will be
                                # reversed compared to what we're looking for
                                if hsp.sbjct_start < hsp.sbjct_end:
                                    end = hsp.sbjct_end
                                else:
                                    end = hsp.sbjct_start
                                # Screen out hits that are shorter than the targets
                                # Keeping this format even though this if statement could be re-written more efficiently
                                if end < alignment.length:
                                    pass
                                else:
                                    # Add the details of the mismatched allele to two dictionaries to be processed below
                                    snpdict[genepath] = hsp.query
                                    datadict[genome] = gene
                            # If alleles aren't being updated, add the allele and the percent identity
                            # match to the reference allele
                            else:
                                self.plusdict[genome][gene][allele] = percentidentity
                        # If the percent identity is below the 98% cutoff threshold or is the hsp is too short
                        elif gene not in self.plusdict[genome] and not snpdict:
                            self.plusdict[genome][gene]["N"] = 0
        # If there are no records, populate the dictionary with "N" and a 0% identity
        else:
            self.plusdict[genome]['{}'.format(gene)]["N"] = 0
        # Add matches that are cutoff < match identity < 99.9% and the same length of the target to the allele file
        if snpdict:
            # Initialise some variables
            allelenames = []  # The allele names already in the allele file
            # Go through the allele files in snpDict
            for gpath in snpdict:
                # Open the allele file
                with open(gpath) as geneFile:
                    for line in geneFile:
                        # Only interested in the header for each allele
                        if ">" in line:
                            # Append all, but only use the last - used to be a string instead of a list
                            allelenames.append(line)
                # Find the allele number and the text before the number for different formats
                allelenumber, alleleprenumber = allelesplitter(allelenames[-1])
                # I wanted to keep the new allele numbers distinct from the old ones, so they will start at 1000000
                if allelenumber < 1000000:
                    newallelenumber = 1000000
                # As this will continuously update the allele database, I need to check if new alleles have been added
                else:
                    newallelenumber = allelenumber + 1
                # Initialise newAllele - formatted updated allele number e.g. >AROC1000000
                newallele = ""
                # Accounts for no pre-number text being present (e.g. >12)
                if not alleleprenumber:
                    # Create a sequence record using BioPython
                    fasta = SeqRecord(Seq(snpdict[gpath], "fasta"),  # snpdict[gpath] is the hsp sequence
                                      description="",  # if this is not added, then some junk is written to the header
                                      id=str(newallelenumber))  # keep the format of the allele number e.g. >1000000
                # If there is pre-number text, do the same thing, but use the allele pre-number as well
                else:
                    alleleprenumber = alleleprenumber.split(">")[1]
                    newallele = '{}_{}'.format(alleleprenumber, newallelenumber)
                    fasta = SeqRecord(Seq(snpdict[gpath], "fasta"),
                                      description="",
                                      id=newallele)
                # Open the allele file to append the new sequence record
                with open(gpath, "a") as updatedFasta:
                    # Use the SeqIO module to properly format the new sequence record
                    SeqIO.write(fasta, updatedFasta, "fasta")
                # Perform some cleanup - need to remove all database and results files associated with the allele file
                # prior to its update. Otherwise, the new changes to the file will not be reflected on subsequent
                # iterations of this genome/allele file comparison
                # Remake the database files
                self.updatedb = [gpath]
                self.makedbthreads(self.updatedb)
                # Now use this new allele in populating self.plusdict
                for updatedGenome in datadict:
                    # The percent identity has to be 100% - this allele matches itself
                    self.plusdict[updatedGenome][datadict[updatedGenome]][newallele] = 100.00

    def sequencetyper(self):
        """Determines the sequence type of each strain based on comparisons to sequence type profiles"""
        # Initialise variables
        header = 0
        # Iterate through the genomes
        for sample in self.metadata:
            genome = sample.name
            # Initialise self.bestmatch[genome] with an int that will eventually be replaced by the number of matches
            self.bestmatch[genome] = defaultdict(int)
            if sample.mlst.profile[self.analysistype] != 'NA':
                # Create the profiledata variable to avoid writing sample.mlst.profiledata[self.analysistype]
                profiledata = sample.mlst.profiledata[self.analysistype]
                # For each gene in plusdict[genome]
                for gene in sample.mlst.genelist[self.analysistype]:
                    # Clear the appropriate count and lists
                    multiallele = []
                    multipercent = []
                    for allele, percentid in self.plusdict[genome][gene].iteritems():
                        # "N" alleles screw up the allele splitter function
                        if allele != "N":
                            # Use the alleleSplitter function to get the allele number
                            allelenumber, alleleprenumber = allelesplitter(allele)
                            # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                            multiallele.append(int(allelenumber))
                            multipercent.append(percentid)
                        # If the allele is "N"
                        else:
                            # Append "N" and a percent identity of 0
                            multiallele.append("N")
                            multipercent.append(0)
                        # Trying to catch cases that where the allele isn't "N", but can't be parsed by alleleSplitter
                        if not multiallele:
                            multiallele.append("N")
                            multipercent.append(0)
                        # Populate self.bestdict with genome, gene, alleles - joined with a space (this was written like
                        # this because allele is a list generated by the .iteritems() above, and the percent identity
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
                            # the result to the bestmatch dictionary. Because genes with multiple alleles were sorted
                            # the same way, these strings with multiple alleles will match: 10 692 will never be 692 10
                            if allele == sortedrefallele:
                                # Increment the number of matches to each profile
                                self.bestmatch[genome][sequencetype] += 1
                # Get the best number of matches
                # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
                try:
                    sortedmatches = sorted(self.bestmatch[genome].items(), key=operator.itemgetter(1), reverse=True)[0][
                        1]
                # If there are no matches, set :sortedmatches to zero
                except IndexError:
                    sortedmatches = 0
                # If there are fewer matches than the total number of genes in the typing scheme
                if 0 < int(sortedmatches) < header:
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
                                # Populate :self.mlstseqtype with the genome, best match to profile, number of matches
                                # to the profile, gene, query allele(s), reference allele(s), and percent identity
                                if self.updateprofile:
                                    self.mlstseqtype[genome][sequencetype][sortedmatches][gene][
                                        str(self.bestdict[genome][gene]
                                            .keys()[0])][sortedrefallele] = str(self.bestdict[genome][gene].values()[0])
                                else:
                                    self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                        self.bestdict[genome][gene]
                                            .keys()[0]] = str(self.bestdict[genome][gene].values()[0])
                    # Add the new profile to the profile file (if the option is enabled)
                    if self.updateprofile:
                        self.reprofiler(int(header), sample.mlst.profile[self.analysistype][0], genome)
                elif sortedmatches == 0:
                    for gene in sample.mlst.genelist[self.analysistype]:
                        # Populate the profile of results with 'negative' values for sequence type and sorted matches
                        self.resultprofile[genome]['0'][sortedmatches][gene][self.bestdict[genome][gene]
                                                                             .keys()[0]] = str(self.bestdict[genome]
                                                                                               [gene].values()[0])
                    # Add the new profile to the profile file (if the option is enabled)
                    if self.updateprofile:
                        self.reprofiler(int(header), sample.mlst.profile[self.analysistype][0], genome)
                # Otherwise, the query profile matches the reference profile
                else:
                    # Iterate through best match
                    for sequencetype, matches in self.bestmatch[genome].iteritems():
                        if matches == sortedmatches:
                            for gene in profiledata[sequencetype]:
                                # Populate resultProfile with the genome, best match to profile, number of matches
                                # to the profile, gene, query allele(s), reference allele(s), and percent identity
                                self.resultprofile[genome][sequencetype][sortedmatches][gene][
                                    self.bestdict[genome][gene]
                                        .keys()[0]] = str(self.bestdict[genome][gene].values()[0])
                dotter()

    def reprofiler(self, header, profilefile, genome):
        # reprofiler(numGenes, profileFile, geneList, genome)
        """
        Creates and appends new profiles as required
        :param header:
        :param profilefile:
        :param genome:
        """
        # Iterate through mlstseqtype - it contains genomes with partial matches to current reference profiles
        # Reset :newprofile
        newprofile = ""
        # Find the last profile entry in the dictionary of profiles
        # Opens uses the command line tool 'tail' to look at the last line of the file (-1). This last line
        # is split on tabs, and only the first entry (the sequence type number) is captured
        profile = subprocess.check_output(['tail', '-1', profilefile]).split("\t")[0]
        # Split the _CFIA from the number - if there is no "_", the just use profile as the profile number
        try:
            profilenumber = int(profile.split("_")[0])
        except IndexError:
            profilenumber = int(profile)
        # If the number is less than 1000000, then new profiles have not previously been added
        if profilenumber < 1000000:
            # Set the new last entry number to be 1000000
            lastentry = 1000000
        # If profiles have previously been added
        else:
            # Set last entry to the highest profile number plus one
            lastentry = profilenumber + 1
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
            newprofile = '{}_CFIA'.format(str(lastentry))
            # The number of matches to the reference profile
            nummatches = self.mlstseqtype[genome][sequencetype].keys()[0]
            for sample in self.metadata:
                if sample.name == genome:
                    # The genes in geneList - should be in the correct order
                    for gene in sample.mlst.genelist[self.analysistype]:
                        # The allele for each gene in the query genome
                        allele = self.mlstseqtype[genome][sequencetype][nummatches][gene].keys()[0]
                        # Append the allele to newprofile
                        newprofile += '\t{}'.format(allele)
                        # Add the MLST results for the query genome as well as the new profile data
                        # to resultProfile
                        self.resultprofile[genome]['{}_CFIA'.format(str(lastentry))][header][gene][allele] = \
                            self.mlstseqtype[genome][sequencetype][nummatches][gene][allele].values()[0]
                    seqcount += 1
        # Only perform the next loop if :newprofile exists
        if newprofile:
            # Open the profile file to append
            with open(profilefile, "ab") as appendfile:
                # Append the new profile to the end of the profile file
                appendfile.write("%s\n" % newprofile)
            # Re-run profiler with the updated files
            self.profiler()

    def reporter(self):
        """ Parse the results into a report"""
        # Initialise variables
        row = ''
        for sample in self.metadata:
            print sample.mlst.reportdir[self.analysistype]
            # make_path(sample.mlst.reportdir[self.analysistype])
            try:

                row += 'Strain,SequenceType,Matches,{},\n'.format(','.join(sample.mlst.genelist[self.analysistype]))
                # print json.dumps(self.resultprofile[sample.name], sort_keys=True, indent=4, separators=(',', ': '))
                # print sample.mlst.
                seqcount = 0
                for seqtype in self.resultprofile[sample.name]:
                    matches = self.resultprofile[sample.name][seqtype].keys()[0]
                    if seqcount == 0:
                        row += '{},{},{},'.format(sample.name, seqtype, matches)
                    else:
                        row += ',{},{},'.format(seqtype, matches)
                    for gene in self.resultprofile[sample.name][seqtype][matches]:
                        allele = self.resultprofile[sample.name][seqtype][matches][gene].keys()[0]
                        percentid = self.resultprofile[sample.name][seqtype][matches][gene].values()[0]
                        # print sample.name, seqtype, matches, gene, allele, percentid
                        row += '{} ({}%),'.format(allele, percentid)
                    row += '\n'
                    seqcount += 1
            except KeyError:
                pass
        # print csvheader
        # print row

    def __init__(self, inputobject):
        self.path = inputobject.path
        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = time.time()
        self.analysistype = inputobject.analysistype
        self.allelefolders = set()
        self.updateallele = inputobject.updateallele
        self.updateprofile = inputobject.updateprofile
        self.updatedb = []
        # Declare queues, list, and dictionaries
        self.dqueue = Queue()
        self.blastqueue = Queue()
        self.plusdict = defaultdict(make_dict)
        self.bestdict = defaultdict(make_dict)
        self.bestmatch = defaultdict(int)
        self.mlstseqtype = defaultdict(make_dict)
        self.resultprofile = defaultdict(make_dict)
        self.profiledata = defaultdict(make_dict)
        # Run the MLST analyses
        self.mlst()
        # print json.dumps(self.plusdict, sort_keys=True, indent=4, separators=(',', ': '))


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


if __name__ == '__main__':
    class Parser(object):

        def strainer(self):
            from accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Not that they must have a file extension that
            # begins with .fa
            self.strains = glob('{}sequences/*.fa*'.format(self.path))
            # Populate the metadata object. This object will be populated to mirror the objects created in the
            # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
            # of a pipeline
            for sample in self.strains:
                # Create the object
                metadata = MetadataObject()
                # Set the base file name of the sequence. Just remove the file extension
                filename = os.path.split(sample)[1].split('.')[0]
                # Set the .name attribute to be the file name
                metadata.name = filename
                # Create the .general attribute
                metadata.general = GenObject()
                #
                metadata.mlst = GenObject()
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = sample
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def organismchooser(self):
            """Allows the user to choose which organism to be used in the analyses"""
            # Initialise a count variable to be used in extracting the desired entry from a list of organisms
            orgcount = 0
            schemecount = 0
            # If the path of the folder containing the allele and profile subfolders is provided
            if self.allelepath:
                # Set the genePath variable for use in blastDatabaseClearer
                self.genepath = '{}alleles'.format(self.allelepath)
                # Remove and previously created blast database files
                blastdatabaseclearer(self.allelepath)
                # Create lists of the alleles, and the profile
                self.alleles = glob('{}/*.*fa*'.format(self.genepath))
                # Get the .txt profile file name and path into a variable
                self.profile = glob('{}profile/*.txt'.format(self.allelepath))
            else:
                # If the name of the organism to analyse was provided
                if not self.organism:
                    # Get a list of the organisms in the (default) Organism subfolder
                    if not self.organismpath:
                        organismlist = glob('{}Organism/*'.format(self.path))
                    elif self.organismpath:
                        organismlist = glob('{}*'.format(self.organismpath))
                    else:
                        organismlist = []
                    # Iterate through the sorted list
                    for folder in sorted(organismlist):
                        # Ensure that folder is, in actuality, a folder
                        if os.path.isdir(folder):
                            # Print out the folder names and the count
                            print "[{}]: {}".format(orgcount, os.path.split(folder)[1])
                            orgcount += 1
                    # Get the user input - the number entered corresponds to the list index
                    response = input("Please select an organism: ")
                    # Get the organism path into a variable
                    organism = sorted(organismlist)[int(response)]
                    self.organism = os.path.split(organism)[1]
                    self.organismpath = self.organismpath if self.organismpath else '{}Organism/{}' \
                        .format(self.path, self.organism)
                # If the name wasn't provided
                else:
                    # Set the organism path as the path + Organism + organism name
                    self.organismpath = '{}Organism/{}'.format(self.path, self.organism)
                if not self.scheme:
                    schemelist = glob('{}/*'.format(self.organismpath))
                    # Iterate through the sorted list
                    for folder in sorted(schemelist):
                        # Ensure that folder is, in actuality, a folder
                        if os.path.isdir(folder):
                            # Print out the folder names and the count
                            print '[{}]: {}'.format(schemecount, os.path.split(folder)[1])
                            schemecount += 1
                    # Same as above
                    schemeresponse = input("Please select a typing scheme:")
                    self.schemepath = sorted(schemelist)[int(schemeresponse)]
                    self.scheme = os.path.split(self.schemepath)[1]
                else:
                    # Otherwise set scheme as follows:
                    self.schemepath = '{}/{}'.format(self.organismpath, self.scheme)
                # Set the variables as above
                self.genepath = '{}/alleles'.format(self.schemepath)
                blastdatabaseclearer(self.genepath)
                # Create lists of the alleles, and the profile
                self.alleles = glob('{}/*.*fa*'.format(self.genepath))
                # Set the name and path of the profile file
                self.profile = glob('{}/profile/*.txt'.format(self.schemepath))
            # Add the appropriate variables to the metadata object for each sample
            for sample in self.samples:
                sample.mlst.alleles = {self.scheme: self.alleles}
                sample.mlst.alleledir = {self.scheme: '{}/'.format(self.genepath)}
                sample.mlst.profile = {self.scheme: self.profile}
                sample.mlst.analysistype = {self.scheme: self.scheme}
                sample.mlst.reportdir = {self.scheme: self.reportpath}

        def __init__(self):
            from argparse import ArgumentParser
            parser = ArgumentParser(description='Performs blast analyses to determine presence of alleles in a genome '
                                                'query, and types genome based on typing profile. Adds novel alleles '
                                                'and profiles to the appropriate files. '
                                                'Example command: '
                                                '-p /home/blais/PycharmProjects/MLST  '
                                                '-s /home/blais/PycharmProjects/MLST/sequences '
                                                '-O /home/blais/PycharmProjects/MLST/Organism '
                                                '-o Vibrio '
                                                '-S MLST')
            parser.add_argument('-p', '--path', required=False, default=os.getcwd(),
                                help='Specify path for custom folder locations. If you don\'t supply additional paths'
                                     'e.g. sequencePath, allelePath, or organismPath, then the program will look for '
                                     'MLST files in .../path/Organism, and the query sequences in ../path/sequences. '
                                     'If you don\'t input a path, then the current working directory will be used.')
            parser.add_argument('-c', '--cutoff', required=False, default=98,
                                help='The percent identity cutoff value for BLAST matches. Default is 98%)')
            parser.add_argument('-s', '--sequencePath', required=False,
                                default='/home/blais/PycharmProjects/MLST/sequences',
                                help='The location of the query sequence files')
            parser.add_argument('-a', '--alleleProfilePath', required=False,
                                # default='/home/blais/PycharmProjects/pythonGeneSeekr/Organism/Salmonella/cgMLST',
                                help='The path of the folder containing the two folders containing '
                                     'the allele files, and the profile file e.g. /folder/path/Organism/Vibrio/cgMLST'
                                     'Please note the requirements for the profile database in the readme')
            parser.add_argument('-O', '--organismPath', required=False,
                                help='The path of the folder containing the organism folders e.g. folder/path/Organism')
            parser.add_argument('-o', '--organism', required=False,
                                help='The name of the organism you wish to type. Must match the folder name containing '
                                     'the schemes e.g. Salmonella')
            parser.add_argument('-S', '--scheme', required=False,
                                help='The scheme you wish to use. Must match the folder name containing the scheme e.g.'
                                     ' cgMLST. Furthermore, this folder must contain two folders: "alleles" and '
                                     '"profile". The alleles folder contains the allele files in .fasta format, and the'
                                     ' profile folder contains the profile in .txt format. Please note the requirements'
                                     ' for the profile in the readme')
            parser.add_argument('-u', '--updateProfileFalse', required=False, action='store_false', default=True,
                                help='By default, the program automatically creates new sequence profiles and appends '
                                     'these profiles to the profile file. If, instead, you wish to wish to see the '
                                     'closest match of a query genome to known reference profiles, set this to False.')
            parser.add_argument('-U', '--updateAlleleFalse', required=False, action='store_false', default=True,
                                help='By default, the program automatically creates new alleles and appends these '
                                     'alleles to the appropriate file. If, instead, you wish to wish to see the '
                                     'closest match of a query genome to known reference alleles, set this to False.')
            parser.add_argument('-r', '--reportDirectory', default='{}/reports'.format(os.getcwd()),
                                help='Path to store the reports defaults to os.getcwd()/reports')

            # Get the arguments into a dictionary
            args = vars(parser.parse_args())
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args['path'], "")
            self.reportpath = os.path.join(args['reportDirectory'], "")
            self.cutoff = float(args['cutoff']) / 100
            if args['sequencePath']:
                self.sequencepath = os.path.join(args['sequencePath'], "")
            else:
                self.sequencepath = ""
            if args['alleleProfilePath']:
                self.allelepath = os.path.join(args['alleleProfilePath'], "")
            else:
                self.allelepath = ""
            if args['organismPath']:
                self.organismpath = os.path.join(args['organismPath'], "")
            else:
                self.organismpath = ""
            self.scheme = args['scheme']
            self.organism = args['organism']
            self.updateprofile = args['updateProfileFalse']
            self.updateallele = args['updateAlleleFalse']

            # Initialise variables
            self.genepath = ''
            self.alleles = ''
            self.profile = ''
            self.strains = []
            self.samples = []
            self.schemepath = ''
            # Get a list of the sequence files
            self.strainer()
            self.organismchooser()


    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.analysistype = self.runmetadata.scheme
            self.alleles = self.runmetadata.alleles
            self.profile = self.runmetadata.profile
            self.cutoff = self.runmetadata.cutoff
            self.updateallele = self.runmetadata.updateallele
            self.updateprofile = self.runmetadata.updateprofile
            self.reportdir = self.runmetadata.reportpath
            # Run the analyses
            MLST(self)

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        from accessoryFunctions import GenObject
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                sample.mlst = GenObject()
                if self.analysistype == 'rmlst':
                    self.alleles = glob('{}rMLST/alleles/*.tfa'.format(self.pipelinefilepath))
                    self.profile = glob('{}rMLST/profile/*.txt'.format(self.pipelinefilepath))
                # Set the metadata file appropriately
                sample.mlst.alleles = {self.analysistype: self.alleles}
                sample.mlst.alleledir = {self.analysistype: '{}rMLST/alleles/'.format(self.pipelinefilepath)}
                sample.mlst.profile = {self.analysistype: self.profile}
                sample.mlst.analysistype = {self.analysistype: self.analysistype}
                sample.mlst.reportdir = {self.analysistype: '{}/{}'.format(sample.general.outputdirectory,
                                                                           self.analysistype)}
            else:
                # Set the metadata file appropriately
                sample.mlst.alleles = {self.analysistype: 'NA'}
                sample.mlst.profile = {self.analysistype: 'NA'}
                sample.mlst.analysistype = {self.analysistype: 'NA'}
                sample.mlst.reportdir = {self.analysistype: 'NA'}

    def __init__(self, inputobject, analysistype):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.pipelinefilepath = inputobject.pipelinefilepath
        self.alleles = ''
        self.profile = ''
        self.cutoff = 100
        self.updateallele = True
        self.updateprofile = True
        # Get the alleles and profile into the metadata
        self.strainer()
        MLST(self)
