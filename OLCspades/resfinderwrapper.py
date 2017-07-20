#!/usr/bin/env python
from glob import glob
from accessoryFunctions import *
import metadataprinter
import subprocess
__author__ = 'adamkoziol'


class Resfinder(object):

    def sequences(self):
        """
        Create objects to store metadata for each strain
        :return: list of metadata objects
        """
        import shutil
        printtime('Finding sequence files', self.start)
        samples = list()
        filelist = glob('{}*.fa*'.format(self.sequencepath))
        # Iterate through the names of the fastq files
        for fasta in sorted(filelist):
            # Set the name
            metadata = MetadataObject()
            name = os.path.split(fasta)[1].split('.')[0]
            metadata.name = name
            # Set the destination folder
            outputdir = os.path.join(self.sequencepath, name)
            # Make the destination folder
            make_path(outputdir)
            # Get the fastq files specific to the fastqname
            specific = glob('{}{}*{}*'.format(self.sequencepath, name, self.extension))
            # Copy the files to the output folder
            try:
                # Link the .gz files to :self.path/:filename
                map(lambda x: shutil.copyfile(x, '{}/{}'.format(outputdir, os.path.split(x)[1])), specific)
            # Except os errors
            except OSError as exception:
                # If there is an exception other than the file exists, raise it
                if exception.errno != errno.EEXIST:
                    raise
            # Initialise the general and run categories
            metadata.general = GenObject()
            metadata.commands = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            metadata.general.bestassemblyfile = glob('{}/*.fa*'.format(outputdir))[0]
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir
            # Attribute added to be compatible with the metadataprinter
            metadata.general.fastqfiles = list()
            # Create the ResFinder genobject
            setattr(metadata, self.analysistype, GenObject())
            metadata[self.analysistype].outputdirectory = outputdir
            samples.append(metadata)
            self.dotter.dotter()
        return samples

    def runner(self):
        """
        Run the appropriate methods in the appropriate order
        """
        # Run resfinder
        self.resfinder()
        metadataprinter.MetadataPrinter(self)

    def resfinder(self):
        """
        Run resfinder on every sample/antimicrobial pair.
        """
        from threading import Thread
        printtime('Performing {} analyses'.format(self.analysistype), self.start)
        # Create the threads for the analysis
        for _ in range(self.cpus):
            threads = Thread(target=self.resthreads, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.runmetadata.samples:
            for antimicrobial in self.antimicrobials:
                self.queue.put((sample, antimicrobial))
        self.queue.join()
        # Reset the dotter counter to 0
        self.dotter.globalcounter()
        # Clean up docker containers
        subprocess.call('docker rm -v $(docker ps -a -q -f status=exited)', shell=True,
                        stdout=self.devnull, stderr=self.devnull)

    def resthreads(self):
        while True:
            sample, antimicrobial = self.queue.get()
            # Set and create the desired output directory for the sample/antimicrobial pair
            outputdir = '{}/{}/{}'.format(sample[self.analysistype].outputdirectory,
                                          self.analysistype, antimicrobial)
            make_path(outputdir)
            # Create the system call to run resfinder.
            call = 'resfinder.pl -d {} -b {} -i {} -a {} -k {} -l {} -o {}'\
                .format(self.targetpath, self.blastpath, sample.general.bestassemblyfile, antimicrobial,
                        self.threshold, self.minlength, outputdir)
            # Only run the analyses if a results file does not already exist
            if not os.path.isfile('{}/results.txt'.format(outputdir)):
                # Run the resfinder call
                subprocess.call(call, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Dots!
            self.dotter.dotter()
            self.queue.task_done()

    def __init__(self, inputobject, analysistype):
        """
        :param inputobject: object containing variables of interest
        """
        from queue import Queue
        import multiprocessing
        # Initialise variables
        self.start = inputobject.starttime
        # Define variables based on supplied arguments
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        self.targetpath = inputobject.targetpath
        self.reportpath = os.path.join(self.path, 'reports')
        self.pipeline = inputobject.pipeline
        self.threshold = inputobject.threshold
        self.minlength = inputobject.minlength
        self.blastpath = inputobject.blastpath
        make_path(self.reportpath)
        self.antimicrobials = ['aminoglycoside', 'beta-lactam', 'colistin', 'fosfomycin', 'fusidicacid', 'macrolide',
                               'nitroimidazole', 'oxazolidinone', 'phenicol', 'quinolone', 'rifampicin', 'sulphonamide',
                               'tetracycline', 'trimethoprim', 'glycopeptide']
        self.analysistype = analysistype
        self.devnull = open(os.devnull, 'wb')
        self.extension = '.fasta'
        self.cpus = multiprocessing.cpu_count()
        self.queue = Queue(maxsize=self.cpus)
        self.dotter = Dotter()
        self.runmetadata = MetadataObject()
        #
        if self.pipeline:
            self.runmetadata.samples = inputobject.runmetadata.samples
            self.runner()
        else:
            self.runmetadata.samples = self.sequences()
            # Reset the dotter counter to 0
            self.dotter.globalcounter()
            self.runner()

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    from .accessoryFunctions import printtime
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify path')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    parser.add_argument('-t', '--targetpath',
                        required=True,
                        help='Path of target files to process.')
    parser.add_argument('-b', '--blastpath',
                        required=True,
                        help='Path to the location with a bin subfolder containing your legacy blast executables e.g. '
                             'supplying /usr will force resfinder to search in /usr/bin')
    parser.add_argument('-k', '--threshold',
                        default=90.00,
                        help='The threshold for % identity for example "95.00" for 95 %')
    parser.add_argument('-l', '--minlength',
                        default=0.60,
                        help='The minimum length of the overlap ex 0.60 for an overlap of minimum 60 %')
    # Get the arguments into an object
    args = parser.parse_args()

    # Create an object to store arguments
    inputs = MetadataObject()
    inputs.path = os.path.join(args.path, '')
    assert os.path.isdir(inputs.path), u'Supplied path is not a valid directory {0!r:s}'.format(inputs.path)
    inputs.sequencepath = os.path.join(args.sequencepath, '')
    assert os.path.isdir(inputs.sequencepath), u'Sequence path  is not a valid directory {0!r:s}' \
        .format(inputs.sequencepath)
    inputs.targetpath = os.path.join(args.targetpath, '')
    inputs.reportpath = os.path.join(inputs.path, 'reports')
    assert os.path.isdir(inputs.targetpath), u'Target path is not a valid directory {0!r:s}' \
        .format(inputs.targetpath)
    inputs.threshold = args.threshold
    inputs.minlength = args.minlength
    inputs.blastpath = args.blastpath
    inputs.starttime = time()
    inputs.pipeline = False
    # Run the pipeline
    Resfinder(inputs, 'ResFinder')
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time() - inputs.starttime) + '\033[0m')
