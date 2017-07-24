#!/usr/bin/env python
from threading import Thread, Lock
from accessoryFunctions import *

__author__ = 'adamkoziol'

threadlock = Lock()


class Prodigal(object):

    def predictthreads(self):
        printtime('Performing gene predictions', self.start)
        # Create the threads for the analyses
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.predict, args=())
                threads.setDaemon(True)
                threads.start()
        for sample in self.metadata:
            # Create the .prodigal attribute
            sample.prodigal = GenObject()
            if sample.general.bestassemblyfile != 'NA':
                self.predictqueue.put(sample)
        self.predictqueue.join()

    def predict(self):
        from subprocess import call
        while True:
            sample = self.predictqueue.get()
            fnull = open(os.devnull, 'w')  # define /dev/null
            reportdir = '{}/prodigal'.format(sample.general.outputdirectory)
            sample.prodigal.reportdir = '{}/prodigal'.format(sample.general.outputdirectory)
            results = '{}/{}_prodigalresults.sco'.format(reportdir, sample.name)
            sample.prodigal.results = results
            prodigal = 'prodigal -i {} -o {} -f sco -d {}/{}_genes.fa'\
                .format(sample.general.bestassemblyfile, results, reportdir, sample.name)
            sample.commands.prodigal = prodigal
            make_path(reportdir)
            size = 0
            if os.path.isfile(results):
                size = os.stat(results).st_size
            if not os.path.isfile(results) or size == 0:
                call(prodigal, shell=True, stdout=fnull, stderr=fnull)
            self.predictqueue.task_done()

    def prodigalparse(self):
        printtime('Parsing gene predictions', self.start)

        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                sample.prodigal.predictedgenestotal = 0
                sample.prodigal.predictedgenesover3000bp = 0
                sample.prodigal.predictedgenesover1000bp = 0
                sample.prodigal.predictedgenesover500bp = 0
                sample.prodigal.predictedgenesunder500bp = 0
                with open(sample.prodigal.results, 'r') as results:
                    for line in results:
                        if line.startswith('>'):
                            start = int(line.split('_')[1])
                            end = int(line.split('_')[2])
                            length = abs(start - end)
                            sample.prodigal.predictedgenestotal += 1
                            if length > 3000:
                                sample.prodigal.predictedgenesover3000bp += 1
                            elif length > 1000:
                                sample.prodigal.predictedgenesover1000bp += 1
                            elif length > 500:
                                sample.prodigal.predictedgenesover500bp += 1
                            else:
                                sample.prodigal.predictedgenesunder500bp += 1

    def __init__(self, inputobject):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.predictqueue = Queue()
        self.parsequeue = Queue()
        self.predictthreads()
        self.prodigalparse()
