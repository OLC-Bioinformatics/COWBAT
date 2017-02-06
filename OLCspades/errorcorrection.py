#!/usr/bin/env python
from accessoryFunctions import *
import quality
__author__ = 'adamkoziol'


class Correct(object):

    def correct(self):
        from threading import Thread
        # Find the fastq files for each sample
        # Only make as many threads are there are samples with fastq files
        for i in range(len([sample.general for sample in self.metadata if type(sample.general.fastqfiles) is list])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.correcting, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the fastqfiles variable - will store trimmed fastq file names if they exist, and raw fastq
            # file names if trimmed fastq files were not created for whatever reason
            if 'trimmedfastqfiles' in sample.general.datastore:
                if type(sample.general.trimmedfastqfiles) is list:
                    fastqfiles = sorted(sample.general.trimmedfastqfiles)
                elif type(sample.general.fastqfiles) is list:
                    fastqfiles = sorted(sample.general.fastqfiles)
                else:
                    fastqfiles = ''
            else:
                fastqfiles = sorted(sample.general.fastqfiles)
            if fastqfiles:
                sample.general.correctedfolder = '{}/corrected'.format(sample.general.outputdirectory)
                if len(fastqfiles) == 1:
                    sample.commands.errorcorrection = 'spades.py --only-error-correction --s1 {} -o {} -t {}'\
                        .format(fastqfiles[0], sample.general.correctedfolder, self.cpus)
                else:
                    sample.commands.errorcorrection = 'spades.py --only-error-correction --pe1-1 {} --pe1-2 {} -o {} ' \
                                                      '-t {}'.format(fastqfiles[0], fastqfiles[1],
                                                                     sample.general.correctedfolder, self.cpus)
                self.correctqueue.put(sample)
        self.correctqueue.join()

    def correcting(self):
        from subprocess import call
        from glob import glob
        while True:
            sample = self.correctqueue.get()
            if not os.path.isdir(sample.general.correctedfolder):
                call(sample.commands.errorcorrection, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Depending on when along pipeline development, analyses were performed, the trimmed, corrected files
            # could be in a different location. Allow for this
            sample.general.correctedfolder = sample.general.correctedfolder \
                if glob('{}/corrected/*_trimmed*'.format(sample.general.correctedfolder)) \
                else '{}/spades_output'.format(sample.general.outputdirectory)
            # Get the trimmed, corrected fastq files into the object
            sample.general.trimmedcorrectedfastqfiles = \
                sorted(glob('{}/corrected/*_trimmed*'.format(sample.general.correctedfolder)))
            self.correctqueue.task_done()

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.cpus = inputobject.cpus
        self.start = inputobject.starttime
        self.correctqueue = Queue(maxsize=self.cpus)
        self.devnull = open(os.devnull, 'wb')
        printtime('Correcting sequences', self.start)
        self.correct()
