#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject
import olctools.accessoryFunctions.metadataprinter as metadataprinter
from genemethods.assemblypipeline import createobject
from argparse import ArgumentParser
from threading import Thread
from subprocess import call
from csv import DictReader
from queue import Queue
import multiprocessing
from time import time
import subprocess
import logging
import os
__author__ = 'adamkoziol'


class FilterGenome(object):

    def objectprep(self):

        # Only find the data files if a datapath is provided
        if self.datapath:
            self.runmetadata = createobject.ObjectCreation(self)
        else:
            for sample in self.runmetadata.samples:
                sample.general.abundancefile = sample.general.abundance
                sample.general.assignmentfile = sample.general.classification
                sample.general.fastqfiles = [sample.general.combined]
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
        # Load the results in the csv files into dictionaries
        self.taxids()

    def taxids(self):
        for sample in self.runmetadata.samples:
            # Initialise a list to store the taxIDs of interest
            sample.general.taxids = list()
            # Read the abundance file into a dictionary
            abundancedict = DictReader(open(sample.general.abundancefile))
            # Filter abundance to taxIDs with at least self.cutoff% of the total proportion
            for row in abundancedict:
                # The UNKNOWN category doesn't contain a 'Lineage' column, and therefore, subsequent columns are
                # shifted out of proper alignment, and do not contain the appropriate data
                try:
                    if float(row['Proportion_All(%)']) > self.cutoff:
                        sample.general.taxids.append(row['TaxID'], )
                except ValueError:
                    pass
            for taxid in sample.general.taxids:
                # Create the an attribute for each taxID
                setattr(sample, taxid, GenObject())
                sample[taxid].readlist = list()
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)
        # Load the assignment file to memory
        self.loadassignment()

    def loadassignment(self):
        """
        Load the taxonomic assignment for each read
        """
        logging.info('Finding taxonomic assignments')
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.assignmentload, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            self.loadqueue.put(sample)
        self.loadqueue.join()
        # Filter the .fastq files
        self.readlist()

    def assignmentload(self):
        while True:
            sample = self.loadqueue.get()
            # Initialise dictionaries to store results
            sample.general.fastqassignment = dict()
            sample.general.reads = dict()
            # Read the assignment file into a dictionary
            with open(sample.general.assignmentfile, 'r') as assignmentcsv:
                for row in assignmentcsv:
                    # Split on ','
                    data = row.split(',')
                    # Each row contains: Object_ID, Length, Assignment
                    # If the Assignment is one of the taxIDs of interest
                    if data[2].rstrip() in sample.general.taxids:
                        # Add the read name to the list of reads associated with that taxID
                        sample[data[2].rstrip()].readlist.append(data[0])
            self.loadqueue.task_done()

    def readlist(self):
        """
        Sort the reads, and create lists to be used in creating sorted .fastq files
        """
        logging.info('Sorting reads')
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.listread, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            self.listqueue.put(sample)
        self.listqueue.join()
        # Create
        self.fastqfilter()

    def listread(self):
        while True:
            sample = self.listqueue.get()
            # Set and create the path of the sorted fastq files
            sample.general.sortedfastqpath = os.path.join(sample.general.outputdirectory, 'sortedFastq')
            make_path(sample.general.sortedfastqpath)
            # Initialise dictionaries to hold data
            sample.general.fastqlist = dict()
            sample.general.filteredfastq = dict()
            # Iterate through the taxIDs
            for taxid in sample.general.taxids:
                # Set the name of the list to store all the reads associated with the taxID
                sample.general.fastqlist[taxid] = os.path.join(sample.general.sortedfastqpath,
                                                               '{sn}_{taxid}.txt'.format(sn=sample.name,
                                                                                         taxid=taxid))
                # Set the name of the .fastq file that will store the filtered reads
                sample.general.filteredfastq[taxid] = os.path.join(sample.general.sortedfastqpath,
                                                                   '{sn}_{taxid}.fastq.gz'.format(sn=sample.name,
                                                                                                  taxid=taxid))
                # Open the list, and write the list of all reads, one per line
                with open(sample.general.fastqlist[taxid], 'w') as binned:
                    binned.write('\n'.join(set(sample[taxid].readlist)))
            self.listqueue.task_done()

    def fastqfilter(self):
        """
        Filter the reads into separate files based on taxonomic assignment
        """
        logging.info('Creating filtered .fastqfiles')
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.filterfastq, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.runmetadata.samples:
            self.filterqueue.put(sample)
        self.filterqueue.join()
        # Print the metadata to file
        metadataprinter.MetadataPrinter(self)

    def filterfastq(self):
        while True:
            sample = self.filterqueue.get()
            # Iterate through the taxIDs
            for taxid in sample.general.taxids:
                # Set the system call to seqtk to subsequence the fastq file
                sample.general.seqtkcall = 'seqtk subseq {fastq} {list} | gzip > {filtered}' \
                    .format(fastq=sample.general.fastqfiles[0],
                            list=sample.general.fastqlist[taxid],
                            filtered=sample.general.filteredfastq[taxid])
                # Run the system call only if the filtered file does not exist
                if not os.path.isfile(sample.general.filteredfastq[taxid]):
                    call(sample.general.seqtkcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                # Delete the large list stored in the object
                delattr(sample[taxid], "readlist")
            self.filterqueue.task_done()

    def __init__(self, inputobject):
        # Define variables based on supplied arguments
        self.start = inputobject.start
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.datapath = inputobject.datapath
        self.reportpath = inputobject.reportpath
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = inputobject.cpus
        # Set the cutoff to be a percent
        self.cutoff = inputobject.cutoff
        # Initialise a variable to hold the sample objects
        self.runmetadata = inputobject.runmetadata if inputobject.runmetadata else MetadataObject()
        # Initialise queues
        self.loadqueue = Queue()
        self.listqueue = Queue()
        self.filterqueue = Queue()
        self.devnull = open(os.devnull, 'wb')


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    class Parser(object):

        def __init__(self):

            # Get the current commit of the pipeline from git
            # Extract the path of the current script from the full path + file name
            homepath = os.path.split(os.path.abspath(__file__))[0]
            # Find the commit of the script by running a command to change to the directory containing the script and
            # run a git command to return the short version of the commit hash
            commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                                      shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
            # Parser for arguments
            parser = ArgumentParser(description='Filter reads based on taxonomic assignment')
            parser.add_argument('-v', '--version',
                                action='version', version='%(prog)s commit {}'.format(commit))
            parser.add_argument('path',
                                help='Specify path')
            parser.add_argument('-t', '--threads',
                                help='Number of threads. Default is the number of cpus in the system')
            parser.add_argument('-s', '--sequencepath',
                                required=True,
                                help='Path of .fastq(.gz) files to process.')
            parser.add_argument('-d', '--datapath',
                                required=True,
                                help='Path of .csv files created by CLARK with read ID, length, and assignment.')
            parser.add_argument('-c', '--cutoff',
                                default=0.01,
                                help='Cutoff value for deciding which taxIDs to use when sorting .fastq files. '
                                     'Defaults to 1 percent. Please note that you must use a decimal format: enter 0.05'
                                     ' to get a 5 percent cutoff value')
            parser.add_argument('-x', '--taxids',
                                help='NOT IMPLEMENTED: CSV of desired taxIDs from each sample. ')
            # Get the arguments into an object
            args = parser.parse_args()
            self.start = time()
            # Define variables based on supplied arguments
            self.path = os.path.join(args.path)
            assert os.path.isdir(self.path), 'Supplied path is not a valid directory {path}'.format(path=self.path)
            self.sequencepath = os.path.join(args.sequencepath)
            assert os.path.isdir(self.sequencepath), 'Sequence location supplied is not a valid directory {seq_path}' \
                .format(seq_path=self.sequencepath)
            self.datapath = os.path.join(args.datapath)
            self.reportpath = os.path.join(self.path, 'reports')
            # Use the argument for the number of threads to use, or default to the number of cpus in the system
            self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
            # Set the cutoff to be a percent
            self.cutoff = args.cutoff * 100
            # Run the pipeline
            self.runmetadata = MetadataObject()
            genome = FilterGenome(self)
            genome.objectprep()
            logging.info('Filtering complete')
    # Run the script
    Parser()


class PipelineInit(object):

    def __init__(self, inputobject):
        # Define variables based on supplied arguments
        self.start = inputobject.start
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.datapath = inputobject.datapath
        self.reportpath = inputobject.reportpath
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = inputobject.cpus
        # Set the cutoff to be a percent
        self.cutoff = inputobject.cutoff
        # Initialise a variable to hold the sample objects
        self.runmetadata = inputobject.runmetadata
        # Initialise queues
        self.loadqueue = Queue()
        self.listqueue = Queue()
        self.filterqueue = Queue()
        self.devnull = open(os.devnull, 'wb')
        # Run the pipeline
        genome = FilterGenome(self)
        genome.objectprep()
        logging.info('Filtering complete')
