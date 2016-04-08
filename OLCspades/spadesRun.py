#!/usr/bin/env python
from accessoryFunctions import *
import os
__author__ = 'adamkoziol'


class Spades(object):

    def spades(self):
        from threading import Thread
        # Find the fastq files for each sample
        # Only make as many threads are there are samples with fastq files
        for i in range(len([sample.general for sample in self.metadata if type(sample.general.fastqfiles) is list])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.assemble, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the spades command
            spadescommand = ''
            # Split the string of the provided kmer argument
            kmerlist = self.kmers.split(',')
            # Regenerate the list of kmers to use if the kmer is less than the readlength
            sample.general.kmers = ','.join([kmer for kmer in kmerlist if int(kmer) <= sample.run.forwardlength])
            # Initialise the fastqfiles variable - will store trimmed fastq file names if they exist, and raw fastq
            # file names if trimmed fastq files were not created for whatever reason
            if type(sample.general.trimmedfastqfiles) is list:
                fastqfiles = sorted(sample.general.trimmedfastqfiles)
            elif type(sample.general.fastqfiles) is list:
                    fastqfiles = sorted(sample.general.fastqfiles)
            else:
                fastqfiles = ''
            # Only proceed if fastq files exists
            if fastqfiles:
                # Set the the forward fastq files
                forward = fastqfiles[0]
                # Set the output directory
                sample.general.spadesoutput = '{}/spades_output'.format(sample.general.outputdirectory)
                # If there are two fastq files
                if len(fastqfiles) == 2:
                    # Set the reverse fastq name
                    reverse = fastqfiles[1]
                    # If a previous assembly was partially completed, continue from the most recent checkpoint
                    if os.path.isdir(sample.general.spadesoutput):
                        spadescommand = 'spades.py -k {} --careful --continue --pe1-1 {} --pe1-2 {} -o {} -t {}'\
                                        .format(sample.general.kmers, forward, reverse, sample.general.spadesoutput,
                                                self.threads)
                    else:
                        spadescommand = 'spades.py -k {} --careful --pe1-1 {} --pe1-2 {} -o {} -t {}'\
                                        .format(sample.general.kmers, forward, reverse, sample.general.spadesoutput,
                                                self.threads)
                # Same as above, but use single read settings for spades
                else:
                    if os.path.isdir(sample.general.spadesoutput):
                        spadescommand = 'spades.py -k {} --careful --continue --s1 {} -o {} -t {}'\
                                        .format(sample.general.kmers, forward, sample.general.spadesoutput,
                                                self.threads)
                    else:
                        spadescommand = 'spades.py -k {} --careful --s1 {} -o {} -t {}'\
                                        .format(sample.general.kmers, forward, sample.general.spadesoutput,
                                                self.threads)
            # If there are no fastq files, populate the metadata appropriately
            else:
                sample.general.spadesoutput = 'NA'
            # Put the arguments to pass to the assemble method into the queue
            self.assemblequeue.put((spadescommand, sample.general.spadesoutput))
            # Add the command to the metadata
            sample.commands.spadescall = spadescommand
            # Add the version to the metadata
            sample.software.spades = get_version(['spades.py']).split('\n')[0].split()[3]
        # Join the threads
        self.assemblequeue.join()
        # Filter contigs shorter than 1000 bp, and rename remaining contigs with sample.name
        self.filter()
        self.insertsize()

    def assemble(self):
        """Run the assembly command in a multi-threaded fashion"""
        from subprocess import call
        while True:
            (command, output) = self.assemblequeue.get()
            if command and not os.path.isfile('{}/contigs.fasta'.format(output)):
                # execute(command)
                call(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            dotter()
            # Signal to the queue that the job is done
            self.assemblequeue.task_done()

    def filter(self):
        """Filter contigs greater than 1000 bp in length, and copy the filtered files to a common assemblies folder"""
        from accessoryFunctions import make_path
        from Bio import SeqIO
        import shutil
        for sample in self.metadata:
            # Set the name of the unfiltered spades assembly output file
            contigsfile = '{}/contigs.fasta'.format(sample.general.spadesoutput)
            # Set the name of the filtered assembly file
            filteredfile = '{}/{}.fasta'.format(sample.general.outputdirectory, sample.name)
            # Only run on samples that have been processed with spades
            if os.path.isfile(contigsfile) and not os.path.isfile(filteredfile):
                # http://biopython.org/wiki/SeqIO#Input.2FOutput_Example_-_Filtering_by_sequence_length
                over1000bp = []
                for record in SeqIO.parse(open(contigsfile, "rU"), "fasta"):
                    # Include only contigs greater than 1000 bp in length
                    if len(record.seq) >= 1000:
                        # Replace 'NODE' in the fasta header with the sample name
                        # >NODE_1_length_705814_cov_37.107_ID_4231
                        # newid = re.sub("NODE", sample.name, record.id)
                        record.id = record.id.replace('NODE', sample.name)
                        # record.id = newid
                        # Clear the name and description attributes of the record
                        record.name = ''
                        record.description = ''
                        # Add this record to our list
                        over1000bp.append(record)
                # Open the filtered assembly file
                with open(filteredfile, 'wb') as formatted:
                    # Write the records in the list to the file
                    SeqIO.write(over1000bp, formatted, 'fasta')
            # If the filtered file was successfully created, copy it to the BestAssemblies folder
            if os.path.isfile(filteredfile):
                # Add the name and path of the filtered file to the metadata
                sample.general.filteredfile = filteredfile
                # Set the assemblies path
                sample.general.bestassembliespath = '{}BestAssemblies'.format(self.path)
                # Make the path (if necessary)
                make_path(sample.general.bestassembliespath)
                # Set the name of the file in the best assemblies folder
                bestassemblyfile = '{}/{}.fasta'.format(sample.general.bestassembliespath, sample.name)
                # Add the name and path of the best assembly file to the metadata
                sample.general.bestassemblyfile = bestassemblyfile
                # Copy the filtered file to the BestAssemblies folder
                if not os.path.isfile(bestassemblyfile):
                    shutil.copyfile(filteredfile, bestassemblyfile)
            else:
                sample.general.bestassemblyfile = 'NA'

    def insertsize(self):
        """Extracts the insert size and its deviation from the spades.log file"""
        for sample in self.metadata:
            # Only look if the spades output folder exists, and if there are two fastq files (can't find the insert
            # size of single reads
            if os.path.isdir(sample.general.spadesoutput) and len(sample.general.fastqfiles) == 2:
                # Set the name of the log file
                spadeslogfile = '{}/spades.log'.format(sample.general.spadesoutput)
                # Open the log file
                with open(spadeslogfile, 'rb') as spadeslog:
                    # Iterate through the file
                    for line in spadeslog:
                        # Find the line with the insert size on it. Will look something like this:
                        """
                        0:02:07.605   144M / 9G    INFO    General (pair_info_count.cpp : 191) \
                        Insert size = 240.514, deviation = 105.257, left quantile = 142, right quantile = 384, \
                        read length = 301
                        """
                        if 'Insert size =' in line:
                            # Extract the relevant data and add it to the metadata
                            sample.general.insertsize = line.split('= ')[1].split(',')[0]
                            sample.general.insertsizestandarddev = line.split('= ')[2].split(',')[0]
            # Otherwise, populate with NA
            else:
                sample.general.insertsize = 'NA'
                sample.general.insertsizestandarddev = 'NA'

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.kmers = inputobject.kmers
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.assemblequeue = Queue()
        printtime('Assembling sequences', self.start)
        self.spades()
