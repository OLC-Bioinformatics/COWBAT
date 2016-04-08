#!/usr/bin/env python
import bz2
from Queue import Queue
from threading import Thread
from glob import glob
from accessoryFunctions import *
__author__ = 'adamkoziol'


class Compress(object):

    def compressthreads(self):
        printtime('Compressing large files', self. start)
        compressfile = []
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                #
                compressfile.append([item for flatlist in
                                     [sample.general.trimmedfastqfiles + sample.general.trimmedcorrectedfastqfiles +
                                      [sample.general.sortedbam]] for item in flatlist])
        #
        compressfiles = [item for flatlist in compressfile for item in flatlist]
        for i in range(len(compressfiles)):
            # Send the threads to makeblastdb
            threads = Thread(target=self.compress, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
            # Make blast databases for MLST files (if necessary)
        for compress in compressfiles:
            self.compressqueue.put(compress)
        self.compressqueue.join()  # wait on the dqueue until everything has been processed
        self.remove()

    def compress(self):
        while True:
            compressfile = self.compressqueue.get()
            if '.bz2' not in compressfile:
                if not os.path.isfile('{}.bz2'.format(compressfile)):
                    output = bz2.BZ2File('{}.bz2'.format(compressfile), 'wb')
                    with open(compressfile, 'rb') as inputfile:
                        for inputdata in inputfile:
                            output.write(inputdata)
                    output.close()
                if os.path.isfile('{}.bz2'.format(compressfile)):
                    try:
                        os.remove(compressfile)
                    except OSError:
                        pass
            self.compressqueue.task_done()

    def remove(self):
        import shutil
        printtime('Removing temporary files', self.start)
        removefile = []
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                removefile.append(glob('{}/K*/'.format(sample.general.spadesoutput)))
        removefiles = [item for flatlist in removefile for item in flatlist]
        for folder in removefiles:
            shutil.rmtree(folder)

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.compressqueue = Queue()
        self.compressthreads()
