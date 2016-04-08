#!/usr/bin/env python
import os
from threading import Thread
from Queue import Queue
import time
from subprocess import call
from accessoryFunctions import printtime
__author__ = 'adamkoziol'


class Quality(object):

    def fastqcthreader(self, level):
        from accessoryFunctions import GenObject, get_version
        from glob import glob
        printtime('Running quality control on {} fastq files'.format(level), self.start)
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                # Send the threads to bbduker. :args is empty as I'm using
                threads = Thread(target=self.fastqc, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            # Create the .software attribute for the metadata
            sample.software = GenObject()
            sample.software.fastqc = get_version(['fastqc', '-v']).split('\n')[0].split()[1]
            fastqccall = ""
            # Check to see if the fastq files exist
            if level == 'Trimmed':
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            elif level == 'trimmedcorrected':
                # Add the location of the corrected fastq files
                sample.general.trimmedcorrectedfastqfiles = sorted(glob('{}/corrected/*trimmed*.gz*'
                                                                        .format(sample.general.spadesoutput)))
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            else:
                fastqfiles = sample.general.fastqfiles
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(fastqfiles) is list:
                # Set the output directory location
                outdir = '{}/fastqc/fastqc{}'.format(sample.general.outputdirectory, level)
                # Separate system calls for paired and unpaired fastq files
                if len(fastqfiles) == 2:
                    # Call fastqc with -q (quiet), -o (output directory), -d (where to store temp files) flags, and
                    # -t (number of threads) flags
                    fastqccall = "fastqc {} {} -q -o {} -t 12".format(fastqfiles[0], fastqfiles[1], outdir)
                elif len(fastqfiles) == 1:
                    fastqccall = "fastqc {} -q -o {} -t 12".format(fastqfiles[0], outdir)
                # Add the arguments to the queue
                sample.commands.fastqccall = fastqccall
                self.qcqueue.put((fastqccall, outdir))
        # Wait on the trimqueue until everything has been processed
        self.qcqueue.join()
        self.qcqueue = Queue()

    def fastqc(self):
        """Run fastqc system calls"""
        from accessoryFunctions import make_path
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, outputdir) = self.qcqueue.get()
            # Check to see if the output directory already exists
            if not os.path.isdir(outputdir):
                # Make the output directory
                make_path(outputdir)
                # Run the system call
                call(systemcall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
                # call(systemcall, shell=True, stdout=devnull, stderr=devnull)
            # Signal to qcqueue that job is done
            self.qcqueue.task_done()

    def trimquality(self):
        """Uses bbduk from the bbmap tool suite to quality and adapter trim"""
        from glob import glob
        print "\r[{:}] Trimming fastq files".format(time.strftime("%H:%M:%S"))
        # Create and start threads for each strain with fastq files
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                # Send the threads to bbduker. :args is empty as I'm using
                threads = Thread(target=self.bbduker, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(sample.general.fastqfiles) is list:
                # Check to see if the fastq files exist
                fastqfiles = sorted(sample.general.fastqfiles)
                # Define the output directory
                outputdir = sample.general.outputdirectory
                # Define the name of the forward trimmed fastq file
                cleanforward = '{}/{}_R1_trimmed.fastq'.format(outputdir, sample.name)
                # Separate system calls for paired and unpaired fastq files
                # TODO minlen=number - incorporate read length
                # http://seqanswers.com/forums/showthread.php?t=42776
                if len(fastqfiles) == 2:
                    cleanreverse = '{}/{}_R2_trimmed.fastq'.format(outputdir, sample.name)
                    bbdukcall = "bbduk.sh -Xmx1g in1={} in2={} out1={} out2={} qtrim=w trimq=20 ktrim=l " \
                        "k=25 mink=11 minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1 tpe tbo" \
                        .format(fastqfiles[0], fastqfiles[1], cleanforward, cleanreverse, self.bbduklocation)
                elif len(fastqfiles) == 1:
                    bbdukcall = "bbduk.sh -Xmx1g in={} out={} qtrim=w trimq=20 ktrim=l k=25 mink=11 " \
                        "minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1" \
                        .format(fastqfiles[0], cleanforward, self.bbduklocation)
                else:
                    bbdukcall = ""
                sample.commands.bbduk = bbdukcall
                # Add the arguments to the queue
                self.trimqueue.put((bbdukcall, cleanforward))
        # Wait on the trimqueue until everything has been processed
        self.trimqueue.join()
        # Add all the trimmed files to the metadata
        for sample in self.metadata:
            # Define the output directory
            outputdir = sample.general.outputdirectory
            # Add the trimmed fastq files to a list
            trimmedfastqfiles = glob('{}/*trimmed.fastq'.format(outputdir, sample.name))
            if not trimmedfastqfiles:
                trimmedfastqfiles = glob('{}/*trimmed.fastq.bz2'.format(outputdir, sample.name))
            # Populate the metadata if the files exist
            sample.general.trimmedfastqfiles = trimmedfastqfiles if trimmedfastqfiles else 'NA'
        print "\r[{:}] Fastq files trimmed".format(time.strftime("%H:%M:%S"))
        self.fastqcthreader('Trimmed')

    def bbduker(self):
        """Run bbduk system calls"""
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, forwardname) = self.trimqueue.get()
            # Check to see if the forward file already exists
            if not os.path.isfile(forwardname) and not os.path.isfile('{}.bz2'.format(forwardname)):
                call(systemcall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            # Signal to trimqueue that job is done
            self.trimqueue.task_done()

    def __init__(self, inputobject):
        from subprocess import Popen, PIPE
        self.metadata = inputobject.runmetadata.samples
        self.qcqueue = Queue()
        self.trimqueue = Queue()
        self.correctqueue = Queue()
        self.start = inputobject.starttime
        # Find the location of the bbduk.sh script. This will be used in finding the adapter file
        self.bbduklocation = os.path.split(Popen('which bbduk.sh', shell=True, stdout=PIPE)
                                           .communicate()[0].rstrip())[0]
