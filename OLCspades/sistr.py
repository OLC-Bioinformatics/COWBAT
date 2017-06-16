#!/usr/bin/env python
from accessoryFunctions import *
__author__ = 'adamkoziol'


class Sistr(object):

    def sistr(self):
        """Perform sistr analyses on Salmonella"""
        from threading import Thread
        printtime('Performing sistr analyses', self.start)
        # Create and start threads
        for i in range(self.cpus):
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.sistrthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            try:
                # Only process strains that have been determined to be Salmonella
                if sample.general.referencegenus == 'Salmonella':
                    # Create the analysis-type specific attribute
                    setattr(sample, self.analysistype, GenObject())
                    # Set and create the path of the directory to store the strain-specific reports
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                    make_path(sample[self.analysistype].reportdir)
                    # Name of the .json output file
                    sample[self.analysistype].jsonoutput = '{}{}.json'\
                        .format(sample[self.analysistype].reportdir, sample.name)
                    # Set the sistr system call
                    sample.commands.sistr = \
                        'sistr -f json -o {} -t 4 -T {}tmp {}'.format(sample[self.analysistype].jsonoutput,
                                                                      sample[self.analysistype].reportdir,
                                                                      sample.general.bestassemblyfile)
                    # Add the sample to the queue
                    self.queue.put(sample)
            except (ValueError, KeyError):
                pass
        self.queue.join()
        self.report()

    def sistrthreads(self):
        from subprocess import call
        import json
        while True:
            # Get the sample object from the queue
            sample = self.queue.get()
            # Only run the analyses if the output json file does not exist
            if not os.path.isfile(sample[self.analysistype].jsonoutput):
                call(sample.commands.sistr, shell=True, stdout=self.devnull, stderr=self.devnull)
            # Read in the output .json file into the metadata
            sample[self.analysistype].jsondata = json.load(open(sample[self.analysistype].jsonoutput, 'rb'))
            self.queue.task_done()

    def report(self):
        """Creates sistr reports"""
        # Initialise strings to store report data
        header = '\t'.join(self.headers) + '\n'
        data = ''
        for sample in self.metadata:
            # Each strain is a fresh row
            row = ''
            try:
                # Set the name of the report.
                # Note that this is a tab-separated file, as there can be commas in the results
                sample[self.analysistype].report = sample[self.analysistype].reportdir + sample.name + '.tsv'
                # Iterate through all the headers to use as keys in the json-formatted output
                for category in self.headers:
                    # Tab separate all the results
                    row += '{}\t'.format(sample[self.analysistype].jsondata[0][category])
                # End the results with a newline
                row += '\n'

                data += row
                # Create and write headers and results to the strain-specific report
                with open(sample[self.analysistype].report, 'wb') as strainreport:
                    strainreport.write(header)
                    strainreport.write(row)
            except (KeyError, AttributeError):
                pass
        # Create and write headers and cumulative results to the combined report
        with open('{}sistr.tsv'.format(self.reportdir), 'wb') as report:
            report.write(header)
            report.write(data)

    def __init__(self, inputobject, analysistype):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.reportdir = '{}/'.format(inputobject.reportpath)
        make_path(self.reportdir)
        self.analysistype = analysistype
        self.devnull = open(os.devnull, 'wb')
        self.queue = Queue()
        self.headers = ['genome', 'cgmlst_distance', 'cgmlst_genome_match', 'cgmlst_matching_alleles', 'h1', 'h2',
                        'serogroup', 'serovar', 'serovar_antigen', 'serovar_cgmlst']
        # Run the analyses
        self.sistr()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    import multiprocessing
    from glob import glob
    # Parser for arguments
    parser = ArgumentParser(description='Automate sistr analyses on a folder of .fasta files')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of .fastq(.gz) files to process.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.starttime = time.time()

    # Find the files
    fastas = sorted(glob(os.path.join(arguments.sequencepath, '*.fa*')))

    # Create a metadata object
    arguments.runmetadata = MetadataObject()
    arguments.runmetadata.samples = list()
    for fasta in fastas:
        metadata = MetadataObject()
        metadata.name = os.path.split(fasta)[1].split('.')[0]
        # Initialise the general and run categories
        metadata.general = GenObject()
        metadata.run = GenObject()
        # Set the destination folder
        outputdir = os.path.join(arguments.sequencepath, metadata.name)
        make_path(outputdir)
        # Add the output directory to the metadata
        metadata.general.outputdirectory = outputdir
        metadata.run.outputdirectory = outputdir
        metadata.general.bestassemblyfile = True
        # Initialise an attribute to store commands
        metadata.commands = GenObject()
        # Assume that all samples are Salmonella
        metadata.general.referencegenus = 'Salmonella'
        # Set the .fasta file as the best assembly
        metadata.general.bestassemblyfile = fasta
        arguments.runmetadata.samples.append(metadata)

    arguments.cpus = multiprocessing.cpu_count()
    arguments.reportpath = os.path.join(arguments.path, 'reports')

    # Run the script
    Sistr(arguments, 'sistr')

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.starttime) + '\033[0m')
