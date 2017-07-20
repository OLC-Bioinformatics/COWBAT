#!/usr/bin/env python
from accessoryFunctions import *
import os
import shutil

__author__ = 'mikeknowles,akoziol'


class Busco(object):
    def buscoprocess(self):
        from threading import Thread
        os.chdir(self.path)
        # Find the fasta files for each sample
        # Only make as many threads are there are samples with fasta files
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.analyze, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Save augustus, blast and BUSCO versions
            sample.software.BUSCO, sample.software.Blastn, sample.software.Augustus, sample.software.python3 = \
                self.version, self.blast, self.augustus, self.pyversion
            if sample.general.bestassemblyfile:
                sample.general.buscoresults = '{}/busco_results'.format(sample.general.outputdirectory)
                buscotemp = "{}run_{}".format(self.path, sample.name)
                sample.commands.BUSCO = "python3 {} -in {} -o {} -l /accessoryfiles/{} -m genome". \
                    format(self.executable, sample.general.bestassemblyfile, sample.name, self.lineage)
                self.qqueue.put((sample, buscotemp))
            else:
                sample.commands.BUSCO = "NA"
        self.qqueue.join()

    def analyze(self):
        """Run the quast command in a multi-threaded fashion"""
        while True:
            sample, temp = self.qqueue.get()
            summary = 'short_summary_{}'.format(sample.name)
            tempfile, moved = [os.path.join(x, summary) for x in [temp, sample.general.buscoresults]]
            # Make sure assembled data exists and BUSCO results do not exist
            if sample.general.bestassemblyfile != 'NA' and map(os.path.isfile, [tempfile, moved]) == [False] * 2:
                if os.path.isdir(temp):  # force incomplete BUSCO runs
                    sample.commands.BUSCO += " -f"
                execute(sample.commands.BUSCO)
            if os.path.isfile(tempfile):
                shutil.move(temp, sample.general.buscoresults)
            if os.path.isfile(moved):
                self.metaparse(sample, moved)
            # Signal to the queue that the job is done
            self.qqueue.task_done()

    @staticmethod
    def metaparse(sample, resfile):
        pc = lambda x: x if x[0].isupper() else x.title()
        if not os.path.isfile(resfile):
            print("There was an issue getting the metadata from {0:s}".format(sample.name))
        else:
            busco = dict()
            # Open BUSCO short_summary file and make list of key value pairs then add those the assembly metadata
            with open(resfile) as report:
                for line in report:
                    # neccesary to split up ifs to avoid exceptions IndexError
                    if line.strip():
                        if line.strip()[0].isdigit():
                            v, k = [[n, "".join([pc(y) for y in k.split()])] for n, k in [line.strip().split('\t')]][0]
                            busco[k] = v
            # TODO: Add support for list of missed BUSCOs
            # This should probably update the datasore to include new busco keyvalue pairs
            sample.assembly.datastore.update(busco)
            # busco.update(sample.assembly.datastore)
            # sample.assembly = GenObject(busco)

    def __init__(self, inputobject):
        from queue import Queue
        from Bio.Blast.Applications import NcbiblastnCommandline
        from distutils import spawn  # TODO: Figure out how to get this imported properly
        # Find blastn and augustus version
        self.version = "v1.1b1"
        self.augustus = " ".join(get_version(['augustus', '--version']).split()[:2])
        self.blast = NcbiblastnCommandline(version=True)()[0].replace('\n', ' ').rstrip()
        self.metadata = inputobject.runmetadata.samples
        # Retrieve abspath of BUSCO executable using spawn
        self.executable = os.path.abspath(spawn.find_executable("BUSCO_{}.py".format(self.version)))
        self.pyversion = get_version(['python3', '-c', 'import sys; print(sys.version)']).rstrip()
        self.start = inputobject.starttime
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.qqueue = Queue()
        printtime('Running BUSCO {} for gene discovery metrics'.format(self.version.split(",")[0]), self.start)
        # Testing with bacterial HMMs
        self.lineage = inputobject.clade
        self.buscoprocess()