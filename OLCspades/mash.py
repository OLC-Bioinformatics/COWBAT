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
                sample.commands.sketch = 'mash sketch -m {} -p {} -l {} -o {}' \
                    .format(self.copies,
                            self.threads,
                            sample[self.analysistype].filelist,
                            sample[self.analysistype].sketchfilenoext)
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
                    'mash dist -p {} {} {} | sort -gk3 > {}'.format(self.threads,
                                                                    sample[self.analysistype].refseqsketch,
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
                # , stdout=self.fnull, stderr=self.fnull
                call(sample.commands.mash, shell=True)
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
        # Allow for the number of threads used to perform the mash sketch and mash dist to decrease as the number
        # of samples increases
        self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        # Use the extension of the files to determine the minimum copies of each k-mer required to pass noise filter
        # in mash sketch; FASTA files: 1 (most of the genome is a single copy, and unique), FASTQ files: 2 (removes
        # singletons, which helps increase accuracy by reducing potential sequencing errors)
        try:
            self.filetype = inputobject.filetype
            self.copies = 2 if 'fastq' in self.filetype or 'gz' in self.filetype else 1
        # Default to FASTQ (2 copies) if the filetype is not provided
        except AttributeError:
            self.copies = 2
        self.sketchqueue = Queue()
        self.mashqueue = Queue()
        self.analysistype = analysistype
        self.fnull = open(os.devnull, 'w')  # define /dev/null
        self.sketching()

if __name__ == '__main__':
    def metadata(sequencepath):
        """
        Creates and populates metadata objects to be consistent with the desired inputs of the Mash class
        :param sequencepath: the path of the folder containing sequence files
        :return: :metadatalist, a list of metadata objects, :extension, the extension of the sequence files
        """
        from accessoryFunctions import MetadataObject, GenObject, make_path, filer
        from glob import glob
        # Get all the sequence files into a list
        files = glob(os.path.join(sequencepath, '*.fa*'))
        # Determine the file extension of the files - this only finds the first extension in the list, so if there are
        # multiple extensions, this will not work
        extension = [os.path.splitext(seqfile)[1] for seqfile in files][0].replace('.', '')
        # Use the filer method to extract the file name from the extension, and any possible additions to the name
        #  e.g. 2013-SEQ-009_S9_L001_R2_001.fastq.gz becomes 2013-SEQ-009
        files = filer(files, extension)
        # Initialise a list to store the metadata objects
        metadatalist = list()
        for seqfile in files:
            # Create the metadata object
            sample = MetadataObject()
            # Populate the metadata object with the required attributes
            sample.name = os.path.split(seqfile)[1]
            sample.general = GenObject()
            sample.commands = GenObject()
            sample.general.bestassemblyfile = glob(os.path.join(sequencepath, sample.name + '*'))
            sample.general.trimmedcorrectedfastqfiles = [seq for seq in sample.general.bestassemblyfile
                                                         if os.path.isfile(seq)]
            sample.general.bestassemblyfile = sample.general.bestassemblyfile[0] if\
                len(sample.general.bestassemblyfile) >= 1 else sample.general.bestassemblyfile
            sample.general.outputdirectory = os.path.join(sequencepath, sample.name)
            # Create the output directory
            make_path(sample.general.outputdirectory)
            # Add the object to the list of objects
            metadatalist.append(sample)
        return metadatalist, extension

    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    import multiprocessing
    # Parser for arguments
    parser = ArgumentParser(description='Run MASH on a directory containing .fasta or .fastq(.gz) files')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of sequences to process with MASH')
    parser.add_argument('-r', '--reffilepath',
                        required=True,
                        help='Path to MASH database')
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Define the start time
    arguments.starttime = time.time()
    arguments.cpus = multiprocessing.cpu_count()
    arguments.reportpath = os.path.join(arguments.path, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create a list of metadata objects, and determine the file extension of the files with the metadata function
    arguments.runmetadata.samples, arguments.filetype = metadata(arguments.sequencepath)
    # Run it
    Mash(arguments, 'mash')
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.starttime) + '\033[0m')
