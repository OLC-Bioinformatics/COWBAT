#!/usr/bin/env python
import os
from glob import glob
from subprocess import call
from threading import Thread
from accessoryFunctions import dotter, printtime
__author__ = 'adamkoziol'


class QualiMap(object):

    def smaltindex(self):
        # Run the indexing threads
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Index the assembly files
            self.indexqueue.put(sample)
        self.indexqueue.join()

    def index(self):
        while True:
            sample = self.indexqueue.get()
            # Set the name of the indexed file name
            filenoext = sample.general.bestassemblyfile.split('.')[0]
            sample.general.filenoext = filenoext
            smifile = filenoext + '.smi'
            # Index the appropriate files if they do not exist
            if not os.path.isfile(smifile):
                # Define the indexing command
                indexcommand = 'cd {} && smalt index {} {}'\
                    .format(sample.general.bestassembliespath, filenoext, sample.general.bestassemblyfile)
                # Run the command
                call(indexcommand, shell=True, stdout=self.fnull, stderr=self.fnull)
            # Print a dot for each indexed file
            dotter()
            # Signal that the thread's task is complete
            self.indexqueue.task_done()

    def smaltmap(self):
        # Run the indexing threads
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.map, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Index the assembly files
            self.mapqueue.put(sample)
        self.mapqueue.join()

    def map(self):
        while True:
            sample = self.mapqueue.get()
            bamfile = sample.general.filenoext + '.bam'
            sample.general.bamfile = bamfile
            # Map the fastq file(s) to the assemblies
            if not os.path.isfile(bamfile):
                # Define the mapping call
                if len(sample.general.trimmedfastqfiles) == 2:
                    smaltmap = 'smalt map -o {} -f bam -n 24 -l pe {} {} {}' \
                               .format(bamfile, sample.general.filenoext, sample.general.trimmedfastqfiles[0],
                                       sample.general.trimmedfastqfiles[1])
                else:
                    smaltmap = 'smalt map -o {} -f bam -n 24 {} {}' \
                               .format(bamfile, sample.general.filenoext, sample.general.trimmedfastqfiles[0])
                # Run the command
                call(smaltmap, shell=True, stdout=self.fnull, stderr=self.fnull)
            # Print a dot for each mapped file
            dotter()
            self.mapqueue.task_done()

    def samtoolssort(self):
        # Run the indexing threads
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.sort, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Index the assembly files
            self.sortqueue.put(sample)
        self.sortqueue.join()

    def sort(self):
        from Bio.Sequencing.Applications import SamtoolsMpileupCommandline
        while True:
            sample = self.sortqueue.get()
            sortedbamfile = sample.general.filenoext + '_sorted.bam'
            sample.general.sortedbamfile = sortedbamfile
            # Map the fastq file(s) to the assemblies
            if not os.path.isfile(sortedbamfile):
                pass
                # Define the sorting call
                # sortcall =
            # Print a dot for each mapped file
            dotter()
            self.sortqueue.task_done()

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.start
        # Define /dev/null
        self.fnull = open(os.devnull, 'wb')
        # Initialise queues
        self.indexqueue = Queue()
        self.mapqueue = Queue()
        self.sortqueue = Queue()
        # Run smalt
        printtime('Indexing assemblies', self.start)
        self.smaltindex()
        printtime('Performing reference mapping', self.start)
        self.smaltmap()
        # Sort the bam files
        printtime('Sorting bam files', self.start)
        self.samtoolssort()


if __name__ == '__main__':
    class Parser(object):

        def associate(self):
            from accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = [fasta for fasta in sorted(glob('{}*.fa*'.format(self.assemblypath)))
                            if '.fastq' not in fasta]
            for strain in self.strains:
                # Extract the name of the strain from the path and file extension
                strainname = os.path.split(strain)[1].split('.')[0]
                # Find the corresponding fastq files for each strain
                fastq = sorted(glob('{}{}*fastq*'.format(self.fastqpath, strainname)))
                # Ensure that fastq files are present for each assembly
                assert fastq, 'Cannot find fastq files for strain {}'.format(strainname)
                # Create the object
                metadata = MetadataObject()
                # Set the .name attribute to be the file name
                metadata.name = strainname
                # Create the .general attribute
                metadata.general = GenObject()
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = strain
                # Set the path of the assembly file
                metadata.general.bestassembliespath = self.assemblypath
                # Populate the .fastqfiles category of :self.metadata
                metadata.general.trimmedfastqfiles = fastq
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            parser = ArgumentParser(description='Calculates coverage depth by mapping FASTQ reads against assemblies')
            parser.add_argument('-p', '--path',
                                default=os.getcwd(),
                                help='Specify the path of the folder that either contains the files of interest, or'
                                     'will be used to store the outputs')
            parser.add_argument('-a', '--assemblies',
                                help='Path to a folder of assemblies. If not provided, the script will look for .fa'
                                     'or .fasta files in the path')
            parser.add_argument('-f', '--fastq',
                                help='Path to a folder of fastq files. If not provided, the script will look for '
                                     'fastq or .fastq.gz files in the path')

            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.assemblypath = os.path.join(args.assemblies, '') if args.assemblies else self.path
            self.fastqpath = os.path.join(args.fastq, '') if args.fastq else self.path

            # Initialise variables
            self.strains = []
            self.samples = []

            # Associate the assemblies and fastq files in a metadata object
            self.associate()

    class MetadataInit(object):
        def __init__(self, start):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.assemblypath = self.runmetadata.assemblypath
            self.fastqpath = self.runmetadata.fastqpath
            self.start = start
            # Run the analyses
            QualiMap(self)

    # Run the class
    from time import time
    starttime = time()
    MetadataInit(starttime)
    printtime('Assembly and characterisation complete', starttime)
