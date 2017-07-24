#!/usr/bin/env python
import os
import time
from queue import Queue
from glob import glob
from subprocess import call
from threading import Thread

from accessoryFunctions import printtime

__author__ = 'adamkoziol'


class Quality(object):

    def fastqcthreader(self, level):
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
            fastqccall = ''
            fastqcreads = ''
            # Check to see if the fastq files exist
            if level == 'Trimmed':

                reader = 'cat' if '.bz2' not in sample.general.trimmedfastqfiles[0] else 'bunzip2 --stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            elif level == 'trimmedcorrected':
                reader = 'gunzip --to-stdout'
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedcorrectedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            else:
                reader = 'cat' if '.gz' not in sample.general.fastqfiles[0] else 'gunzip --to-stdout'
                fastqfiles = sample.general.fastqfiles
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(fastqfiles) is list:
                # Set the output directory location
                outdir = '{}/fastqc/fastqc{}'.format(sample.general.outputdirectory, level)
                # Separate system calls for paired and unpaired fastq files
                if len(fastqfiles) == 2:
                    # Combine the fastq files, so that the paired files are processed together
                    # Call fastqc with -q (quiet), -o (output directory), -t (number of threads) flags
                    fastqccall = '{} {} {} | fastqc -q -t {} stdin -o {}'\
                        .format(reader, fastqfiles[0], fastqfiles[1], self.cpus, outdir)
                    # Also perform QC on the individual forward and reverse reads
                    fastqcreads = "fastqc {} {} -q -o {} -t 12".format(fastqfiles[0], fastqfiles[1], outdir)
                elif len(fastqfiles) == 1:
                    # fastqccall = "fastqc {} -q -o {} -t 12".format(fastqfiles[0], outdir)
                    fastqccall = '{} {} | fastqc -q -t {} stdin -o {}'\
                        .format(reader, fastqfiles[0], self.cpus, outdir)
                    fastqcreads = "fastqc {} -q -o {} -t 12".format(fastqfiles[0], outdir)
                # Add the arguments to the queue
                sample.commands.fastqc = fastqccall
                self.qcqueue.put((sample, fastqccall, outdir, fastqcreads))
        # Wait on the trimqueue until everything has been processed
        self.qcqueue.join()
        self.qcqueue = Queue()

    def fastqc(self):
        """Run fastqc system calls"""
        from accessoryFunctions import make_path
        import shutil
        while True:  # while daemon
            # Unpack the variables from the queue
            (sample, systemcall, outputdir, fastqcreads) = self.qcqueue.get()
            # Check to see if the output HTML file already exists
            if not os.path.isfile('{}/{}_fastqc.html'.format(outputdir, sample.name)):
                # Make the output directory
                make_path(outputdir)
                # Run the system calls
                call(systemcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                call(fastqcreads, shell=True, stdout=self.devnull, stderr=self.devnull)
                # Rename the outputs
                try:
                    shutil.move('{}/stdin_fastqc.html'.format(outputdir),
                                '{}/{}_fastqc.html'.format(outputdir, sample.name))
                    shutil.move('{}/stdin_fastqc.zip'.format(outputdir),
                                '{}/{}_fastqc.zip'.format(outputdir, sample.name))
                except IOError:
                    pass
            # Signal to qcqueue that job is done
            self.qcqueue.task_done()

    def trimquality(self):
        """Uses bbduk from the bbmap tool suite to quality and adapter trim"""
        printtime("Trimming fastq files", self.start)
        # print("\r[{:}] Trimming fastq files".format(time.strftime("%H:%M:%S")))
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
                # Define the name of the trimmed fastq files
                cleanforward = '{}/{}_R1_trimmed.fastq'.format(outputdir, sample.name)
                cleanreverse = '{}/{}_R2_trimmed.fastq'.format(outputdir, sample.name)
                if self.numreads == 2:
                    # Separate system calls for paired and unpaired fastq files
                    # TODO minlen=number - incorporate read length
                    # http://seqanswers.com/forums/showthread.php?t=42776
                    # BBduk 37.23 doesn't need the ktrim=l/mink=11 parameters, so they have been removed.
                    if len(fastqfiles) == 2:
                        if int(sample.run.forwardlength) > 75 and int(sample.run.reverselength) > 75:
                            bbdukcall = "bbduk2.sh -Xmx1g in1={} in2={} out1={} out2={} qtrim=w trimq=20 " \
                                "k=25 minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1 " \
                                        "tpe tbo" \
                                .format(fastqfiles[0], fastqfiles[1], cleanforward, cleanreverse, self.bbduklocation)
                        else:
                            bbdukcall = "bbduk2.sh -Xmx1g in1={} out1={} qtrim=w trimq=20 k=25 " \
                                        "minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1" \
                                .format(fastqfiles[1], cleanreverse, self.bbduklocation)
                    elif len(fastqfiles) == 1:
                        bbdukcall = "bbduk2.sh -Xmx1g in={} out={} qtrim=w trimq=20 k=25 " \
                            "minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1" \
                            .format(fastqfiles[0], cleanforward, self.bbduklocation)
                    else:
                        bbdukcall = ""
                # Allows for exclusion of the reverse reads if desired
                else:
                    bbdukcall = "bbduk2.sh -Xmx1g in={} out={} qtrim=w trimq=20 k=25 " \
                                "minlength=50 forcetrimleft=15 ref={}/resources/adapters.fa hdist=1" \
                        .format(fastqfiles[0], cleanforward, self.bbduklocation)
                    # There is a check to ensure that the trimmed reverse file is created. This will change the file
                    # being looked for to the forward file
                    cleanreverse = cleanforward
                    if self.forwardlength != 'full':
                        bbdukcall += ' forcetrimright={}'.format(str(self.forwardlength))
                sample.commands.bbduk = bbdukcall
                # Add the arguments to the queue
                self.trimqueue.put((sample, bbdukcall, cleanreverse))
        # Wait on the trimqueue until everything has been processed
        self.trimqueue.join()
        # Add all the trimmed files to the metadata

        #print("\r[{:}] Fastq files trimmed".format(time.strftime("%H:%M:%S")))
        printtime('Fastq files trimmed', self.start)
        self.fastqcthreader('Trimmed')

    def bbduker(self):
        """Run bbduk system calls"""
        while True:  # while daemon
            # Unpack the variables from the queue
            (sample, systemcall, reversename) = self.trimqueue.get()
            # Check to see if the forward file already exists
            if systemcall:
                if not os.path.isfile(reversename) and not os.path.isfile('{}.bz2'.format(reversename)):
                    call(systemcall, shell=True, stdout=self.devnull, stderr=self.devnull)
                # Define the output directory
                outputdir = sample.general.outputdirectory
                # Add the trimmed fastq files to a list
                trimmedfastqfiles = sorted(glob('{}/*trimmed.fastq'.format(outputdir, sample.name)))
                if not trimmedfastqfiles:
                    trimmedfastqfiles = sorted(glob('{}/*trimmed.fastq.bz2'.format(outputdir, sample.name)))
                # Populate the metadata if the files exist
                sample.general.trimmedfastqfiles = trimmedfastqfiles if trimmedfastqfiles else 'NA'
            # Signal to trimqueue that job is done
            self.trimqueue.task_done()

    def __init__(self, inputobject):
        from subprocess import Popen, PIPE
        self.metadata = inputobject.runmetadata.samples
        self.cpus = inputobject.cpus
        self.devnull = open(os.devnull, 'wb')
        self.qcqueue = Queue(maxsize=self.cpus)
        self.trimqueue = Queue(maxsize=self.cpus)
        self.correctqueue = Queue(maxsize=self.cpus)
        self.start = inputobject.starttime
        self.forwardlength = inputobject.forwardlength
        self.reverselength = inputobject.reverselength
        self.numreads = inputobject.numreads
        # Find the location of the bbduk.sh script. This will be used in finding the adapter file
        self.bbduklocation = os.path.split(Popen('which bbduk.sh', shell=True, stdout=PIPE)
                                           .communicate()[0].rstrip())[0]
        self.bbduklocation = self.bbduklocation.decode('utf-8')
