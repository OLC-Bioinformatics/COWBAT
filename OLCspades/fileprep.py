#!/usr/bin/env python
import gzip
import os
from threading import Thread
from queue import Queue
__author__ = 'adamkoziol'


class Fileprep(object):

    def fileprep(self):
        """Decompress and concatenate .fastq files"""
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.prep, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Set the name of the decompressed, combined .fastq file
            sample.general.combined = '{}/{}_combined.fastq'.format(sample.general.outputdirectory, sample.name)
            self.queue.put(sample)
        self.queue.join()

    def prep(self):
        while True:
            sample = self.queue.get()
            # Don't make the file if it already exists
            if not os.path.isfile(sample.general.combined):
                # Open this .fastq file to write all the decompressed reads
                #with open(sample.general.combined, 'w') as combined:
                    # Iterate through the uncompressed .fastq file(s)
                for read in sample.general.fastqfiles:
                    # Only decompress if the reads are gzipped
                    if '.gz' in read:
                        with open(sample.general.combined, 'wb') as combined:
                            # Open the .fastq file with gzip
                            with gzip.open(read, 'rb') as fastq:
                                # Read the file contents and write them to the combined file
                                combined.write(fastq.read())
                    else:
                        with open(sample.general.combined, 'w') as combined:
                            with open(read, 'r') as fastq:
                                # Read in data and write it to file
                                combined.write(fastq.read())
            self.queue.task_done()

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.cpus = inputobject.cpus
        self.queue = Queue()
        # Prep the files
        self.fileprep()
