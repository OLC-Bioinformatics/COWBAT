#!/usr/bin/env python
from subprocess import call
from threading import Thread

from accessoryFunctions import *

__author__ = 'adamkoziol'


class Mash(object):
    def sketching(self):
        printtime('Indexing assemblies for mash analysis', self.starttime)
        # Create the threads for the analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.sketch, args=())
                threads.setDaemon(True)
                threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            # Create the analysis type-specific GenObject
            setattr(sample, self.analysistype, GenObject())
            if sample.general.bestassemblyfile != 'NA':
                # Set attributes
                sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory, self.analysistype)
                sample[self.analysistype].targetpath = os.path.join(self.referencefilepath, self.analysistype)
                sample[self.analysistype].refseqsketch = \
                    sample[self.analysistype].targetpath + '/RefSeqSketchesDefaults.msh'
                sample[self.analysistype].sketchfilenoext = '{}/{}'.format(sample[self.analysistype].reportdir,
                                                                           sample.name)
                sample[self.analysistype].sketchfile = sample[self.analysistype].sketchfilenoext + '.msh'
                # Make the mash output directory if necessary
                make_path(sample[self.analysistype].reportdir)
                # Create a file containing the path/name of the filtered, corrected fastq files
                sample[self.analysistype].filelist = '{}/{}_fastqfiles.txt'.format(sample[self.analysistype].reportdir,
                                                                                   sample.name)
                with open(sample[self.analysistype].filelist, 'w') as filelist:
                    filelist.write('\n'.join(sample.general.trimmedcorrectedfastqfiles))

                # Create the system call
                sample.commands.sketch = 'mash sketch -m 2 -p {} -l {} -o {}' \
                    .format(self.cpus, sample[self.analysistype].filelist, sample[self.analysistype].sketchfilenoext)
                # Add each sample to the threads
                self.sketchqueue.put(sample)
        # Join the threads
        self.sketchqueue.join()
        self.mashing()

    def sketch(self):
        while True:
            sample = self.sketchqueue.get()
            if not os.path.isfile(sample[self.analysistype].sketchfile):
                call(sample.commands.sketch, shell=True, stdout=self.fnull, stderr=self.fnull)
            self.sketchqueue.task_done()

    def mashing(self):
        printtime('Performing mash analyses', self.starttime)
        # Create the threads for the analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.mash, args=())
                threads.setDaemon(True)
                threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                sample[self.analysistype].mashresults = '{}/{}.tab'.format(sample[self.analysistype].reportdir,
                                                                           sample.name)

                sample.commands.mash = \
                    'mash dist -p {} {} {} | sort -gk3 > {}'.format(self.cpus, sample[self.analysistype].refseqsketch,
                                                                    sample[self.analysistype].sketchfile,
                                                                    sample[self.analysistype].mashresults)
                self.mashqueue.put(sample)
        # Join the threads
        self.mashqueue.join()
        self.parse()

    def mash(self):
        while True:
            sample = self.mashqueue.get()
            if not os.path.isfile(sample[self.analysistype].mashresults):
                call(sample.commands.mash, shell=True, stdout=self.fnull, stderr=self.fnull)
            self.mashqueue.task_done()

    def parse(self):
        import re
        printtime('Determining closest refseq genome', self.starttime)
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Open the results and extract the first line of data
                mashdata = open(sample[self.analysistype].mashresults).readline().rstrip()
                # Split on tabs
                data = mashdata.split('\t')
                referenceid, queryid, sample[self.analysistype].mashdistance, sample[self.analysistype]. \
                    pvalue, sample[self.analysistype].nummatches = data
                # The database is formatted such that the reference file name is usually preceded by '-.-'
                # e.g. refseq-NZ-1005511-PRJNA224116-SAMN00794588-GCF_000303935.1-.-Escherichia_coli_PA45.fna
                #      refseq-NZ-1639-PRJNA224116-SAMN03349770-GCF_000951975.1-p3KSM-Listeria_monocytogenes.fna
                sample[self.analysistype].closestrefseq = re.findall(r'.+-(.+)\.fna', referenceid)[0]
                # try:
                #     print re.findall(r'.+-(.+)\.fna', referenceid), referenceid
                #     quit()
                # except AttributeError:
                #     print referenceid
                #     quit()
                # try:
                #     sample[self.analysistype].closestrefseq = \
                #         re.search('(?:GCF_.{11}-.-)(.+)\.fna', referenceid).groups()[0]
                # except AttributeError:
                #     try:
                #         sample[self.analysistype].closestrefseq = \
                #             referenceid.split('-.-')[1].split('.fna')[0]
                #     except IndexError:
                #         print sample.name
                #         sample[self.analysistype].closestrefseq = 'NA'
                sample[self.analysistype].closestrefseqgenus = sample[self.analysistype].closestrefseq.split('_')[0]
            else:
                # Populate the attribute with negative results
                sample[self.analysistype].closestrefseqgenus = 'NA'
        # Create the report
        self.reporter()

    def reporter(self):
        make_path(self.reportpath)
        header = 'Strain,ReferenceGenus,ReferenceFile,ReferenceGenomeMashDistance,Pvalue,NumMatchingHashes\n'
        data = ''
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                data += '{},{},{},{},{},{}\n'.format(sample.name,
                                                     sample[self.analysistype].closestrefseqgenus,
                                                     sample[self.analysistype].closestrefseq,
                                                     sample[self.analysistype].mashdistance,
                                                     sample[self.analysistype].pvalue,
                                                     sample[self.analysistype].nummatches)
        # Create the report file
        reportfile = '{}/mash.csv'.format(self.reportpath)
        with open(reportfile, 'w') as report:
            report.write(header)
            report.write(data)

    def __init__(self, inputobject, analysistype):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.referencefilepath = inputobject.reffilepath
        self.starttime = inputobject.starttime
        self.reportpath = inputobject.reportpath
        self.cpus = inputobject.cpus
        self.sketchqueue = Queue()
        self.mashqueue = Queue()
        self.analysistype = analysistype
        self.fnull = open(os.devnull, 'w')  # define /dev/null
        self.sketching()
