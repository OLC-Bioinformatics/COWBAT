#!/usr/bin/env python
from accessoryFunctions import *
__author__ = 'adamkoziol'


class Quast(object):

    def quast(self):
        from threading import Thread
        printtime('Performing Quast analyses', self.start)
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile != 'NA'])):
            # Send the threads to the merge method. :args is empty
            threads = Thread(target=self.runquast, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Create the quast output directory
                quastoutputdirectory = '{}/quast_results/'.format(sample.general.outputdirectory)
                make_path(quastoutputdirectory)
                # If the best reference genome was identified using rMLST, perform the GAGE analysis
                # if sample['rmlst'].referencegenomepath != 'NA':
                #     quastcall = 'quast.py -R {} --gage {} -o {}'.format(sample['rmlst'].referencegenomepath,
                #                                                         sample.general.filteredfile,
                #                                                         quastoutputdirectory)
                # Otherwise run quast without GAGE analyses
                # else:
                quastcall = 'quast.py {} -o {}'.format(sample.general.filteredfile, quastoutputdirectory)
                # Add the command to the metadata
                sample.commands.quast = quastcall
                self.quastqueue.put((sample, quastoutputdirectory))
            else:
                sample.commands.quast = 'NA'
        self.quastqueue.join()

    def runquast(self):
        from subprocess import call
        while True:
            sample, quastoutputdirectory = self.quastqueue.get()
            make_path(quastoutputdirectory)
            fnull = open(os.devnull, 'wb')
            # Don't re-perform the analysis if the report file exists
            if not os.path.isfile('{}/report.tsv'.format(quastoutputdirectory)):
                call(sample.commands.quast, shell=True, stdout=fnull, stderr=fnull)
            # Following the analysis, parse the report (if it exists) into the metadata object
            if os.path.isfile('{}/report.tsv'.format(quastoutputdirectory)):
                self.metaparse(sample, quastoutputdirectory)
            self.quastqueue.task_done()

    def metaparse(self, sample, quastoutputdirectory):
        import functools
        # Tuples of strings to replace when parsing the results file
        repls = ('>=', 'Over'), ('000 Bp', 'kbp'), ('#', 'Num'), \
                ("'", ''), ('(', ''), (')', ''), (' ', ''), ('>', 'Less'), ('Gc%', 'GC%')
        # Initialise the results dictionary
        quast = dict()
        # The results file is gage_report.tsv if that file exists, otherwise it is report.tsv
        resfile = "{0:s}/gage_report.tsv".format(quastoutputdirectory) \
            if os.path.isfile("{0:s}/gage_report.tsv".format(quastoutputdirectory)) \
            else "{0:s}/report.tsv".format(quastoutputdirectory)
        with open(resfile) as report:
            report.next()
            for line in report:
                # Use headings in report as keys for the GenObject supplied from generator and replace incrementally
                # with reduce and lambda function below
                k, v = [functools.reduce(lambda a, kv: a.replace(*kv), repls, s.title()) for s in line.rstrip().split('\t')]
                quast[k] = v
        # Create the quast metadata object
        sample.quast = GenObject(quast)
        sample.quast.outputdirectory = quastoutputdirectory
        sample.quast.kmers = self.kmers

    def __init__(self, inputobject):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.kmers = inputobject.kmers
        self.start = inputobject.starttime
        self.quastqueue = Queue()
        self.quast()
