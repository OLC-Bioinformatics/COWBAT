#!/usr/bin/env python
from glob import glob
from subprocess import call
from threading import Thread
from accessoryFunctions import *
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
__author__ = 'mike knowles, adamkoziol'


class QualiMap(object):

    def smaltindex(self):
        # Run the indexing threads
        for i in range(len([sample.general for sample in self.metadata if sample.general.filteredfile])):
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
            filenoext = sample.general.filteredfile.split('.')[0]
            sample.general.filenoext = filenoext
            smifile = filenoext + '.smi'
            # Define the indexing command
            indexcommand = 'cd {} && smalt index {} {}'\
                .format(sample.general.bestassembliespath, filenoext, sample.general.filteredfile)
            # Index the appropriate files if they do not exist
            if not os.path.isfile(smifile):
                # Run the command
                call(indexcommand, shell=True, stdout=self.fnull, stderr=self.fnull)
            sample.commands.smaltindex = indexcommand
            # Print a dot for each indexed file
            dotter()
            # Signal that the thread's task is complete
            self.indexqueue.task_done()

    def smaltmap(self):
        # Run the indexing threads
        for i in range(len([sample.general for sample in self.metadata if sample.general.filteredfile])):
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
            bamfile = sample.general.filenoext + '_sorted'
            sample.mapping.BamFile = sample.general.filenoext + '_sorted.bam'
            sample.general.bamfile = bamfile
            # Map the fastq file(s) to the assemblies
            # Define the mapping call
            if len(sample.general.fastqfiles) == 2:
                # Paired-end system call. Note that the output from SMALT is piped into samtools sort to prevent
                # the creation of intermediate, unsorted bam files
                smaltmap = 'smalt map -f bam -n {} -x {} {} {} | samtools sort - {}' \
                           .format(self.cpus, sample.general.filenoext, sample.general.fastqfiles[0],
                                   sample.general.fastqfiles[1], bamfile)
            else:
                smaltmap = 'smalt map -f bam -n {} {} {} | samtools sort - {}' \
                           .format(self.cpus, sample.general.filenoext, sample.general.fastqfiles[0],
                                   bamfile)
            # Populate metadata
            sample.software.SMALT = self.smaltversion
            sample.software.samtools = self.samversion
            sample.commands.smaltsamtools = smaltmap
            # Run the call if the sorted bam file doesn't exist
            size = 0
            if os.path.isfile(sample.mapping.BamFile):
                size = os.stat(sample.mapping.BamFile[0]).st_size
            if not os.path.isfile(sample.mapping.BamFile) or size == 0:
                # Run the command
                call(smaltmap, shell=True, stdout=self.fnull, stderr=self.fnull)
            # Print a dot for each mapped file
            dotter()
            self.mapqueue.task_done()

    def __call__(self):
        """Execute Qualimap on call"""
        printtime('Reading BAM file for Qualimap output', self.start)
        for i in range(len([sample.general for sample in self.metadata if sample.general.filteredfile != "NA"])):
            # Send the threads to the merge method. :args is empty
            threads = Thread(target=self.mapper, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.filteredfile != "NA":
                # Set the results folder
                sample.general.QualimapResults = '{}/qualimap_results'.format(sample.general.outputdirectory)
                # Create this results folder if necessary
                make_path(sample.general.QualimapResults)
                sample.software.Qualimap = self.version
                # Define the Qualimap call
                sample.commands.Qualimap = 'qualimap bamqc -bam {} -outdir {}'. \
                    format(sample.mapping.BamFile, sample.general.QualimapResults)
                self.qqueue.put(sample)
            else:
                sample.commands.Qualimap = "NA"
        self.qqueue.join()

    def mapper(self):
        while True:
            sample = self.qqueue.get()
            if sample.general.filteredfile != "NA":
                # Define the Qualimap log and report files
                log = os.path.join(sample.general.QualimapResults, "qualimap.log")
                reportfile = os.path.join(sample.general.QualimapResults, 'genome_results.txt')
                # Initialise a dictionary to hold the Qualimap results
                qdict = dict()
                # If the report file doesn't exist, run Qualimap, and print logs to the log file
                if not os.path.isfile(reportfile):
                    execute(sample.commands.Qualimap, log)
                # Otherwise open the report
                else:
                    with open(reportfile) as report:
                        # Read the report
                        for line in report:
                            # Sanitise the keys and values using self.analyze
                            key, value = self.analyze(line)
                            # If the keys and values exist, enter them into the dictionary
                            if (key, value) != (None, None):
                                qdict[key] = value
                # If there are values in the dictionary
                if qdict:
                    # Make new category for Qualimap results and populate this category with the report data
                    setattr(sample, "mapping", GenObject(qdict))
            self.qqueue.task_done()

    @staticmethod
    def analyze(line):
        # Split on ' = '
        if ' = ' in line:
            key, value = line.split(' = ')
            # Replace occurrences of
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Should we keep comma separation?
            value = value.replace(",", "").replace(" ", "").rstrip()
        # Otherwise set the keys and values to None
        else:
            key, value = None, None
        return key, value

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        # Define /dev/null
        self.fnull = open(os.devnull, 'wb')
        self.smaltversion = get_version(['smalt', 'version']).split('\n')[2].split()[1]
        self.samversion = get_version(['samtools']).split('\n')[2].split()[1]
        self.version = get_version(['qualimap', '--help']).split('\n')[4].split()[1]
        # Initialise queues
        self.indexqueue = Queue()
        self.mapqueue = Queue()
        self.sortqueue = Queue()
        self.qqueue = Queue()
        # Run smalt
        printtime('Indexing assemblies', self.start)
        self.smaltindex()
        printtime('Performing reference mapping', self.start)
        self.smaltmap()


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
                # Set the .general.filteredfile file to be the name and path of the sequence file
                metadata.general.filteredfile = strain
                # Set the path of the assembly file
                metadata.general.bestassembliespath = self.assemblypath
                # Populate the .fastqfiles category of :self.metadata
                metadata.general.trimmedfastqfiles = fastq
                # Create the output directory path
                metadata.general.outputdirectory = '{}{}'.format(self.path, strainname)
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            import subprocess
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
            parser.add_argument('-t', '--threads',
                                help='Number of threads. Default is the number of cores in the system')
            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.assemblypath = os.path.join(args.assemblies, '') if args.assemblies else self.path
            self.fastqpath = os.path.join(args.fastq, '') if args.fastq else self.path
            # Use the argument for the number of threads to use, or default to the number of cpus in the system
            self.cpus = args.threads if args.threads else int(subprocess.Popen("awk '/^processor/ { N++} END "
                                                                               "{ print N }' /proc/cpuinfo",
                                                                               shell=True,
                                                                               stdout=subprocess.PIPE)
                                                              .communicate()[0].rstrip())
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
            self.starttime = start
            self.cpus = self.runmetadata.cpus
            # Run the analyses - the extra set of parentheses is due to using the __call__ method in the class
            QualiMap(self)()

    # Run the class
    from time import time
    starttime = time()
    MetadataInit(starttime)
    printtime('Assembly and characterisation complete', starttime)
