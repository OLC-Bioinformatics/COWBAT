#!/usr/bin/env python
from threading import Thread
from accessoryFunctions import *
import metadataprinter

__author__ = 'adamkoziol'


class Vtyper(object):

    def vtyper(self):
        """Setup and create  threads for ePCR"""
        printtime('Running ePCR', self.start)
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.epcr, args=())
                threads.setDaemon(True)
                threads.start()
        # Create the system calls for famap, fahash, and ePCR
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'stx' in sample.general.datastore:
                    setattr(sample, self.analysistype, GenObject())
                    # Get the primers ready
                    if self.reffilepath:
                        sample[self.analysistype].primers = '{}{}/vtx_subtyping_primers.txt'\
                            .format(self.reffilepath, self.analysistype)
                    else:
                        sample[self.analysistype].primers = self.primerfile
                    # Make the output path
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                    make_path(sample[self.analysistype].reportdir)
                    outfile = sample[self.analysistype].reportdir + sample.name
                    # Set the hashing and mapping commands
                    sample.commands.famap = 'famap -b {}.famap {}.fasta'.format(outfile, sample.general.filenoext)
                    sample.commands.fahash = 'fahash -b {}.hash {}.famap'.format(outfile, outfile)
                    # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
                    # -S {hash file} (Perform STS lookup using hash-file),
                    # -r + (Enable/disable reverse STS lookup)
                    # -m 10000 (Set variability for STS size for lookup),
                    # -n 1 (Set max allowed mismatches per primer for lookup)
                    # -g 0 (Set max allowed indels per primer for lookup),
                    # -G (Print alignments in comments),
                    # -q quiet
                    # -o {output file},
                    sample.commands.epcr = 're-PCR -S {}.hash -r + -m 10000 -n 1 -g 0 -G -q -o {}.txt {}'\
                        .format(outfile, outfile, sample[self.analysistype].primers)
                    sample[self.analysistype].resultsfile = '{}.txt'.format(outfile)
                    self.epcrqueue.put((sample, outfile))
        self.epcrqueue.join()
        self.epcrparse()

    def epcr(self):
        from subprocess import call
        while True:
            sample, linkfile = self.epcrqueue.get()
            if not os.path.isfile('{}.famap'.format(linkfile)):
                call(sample.commands.famap, shell=True, stdout=self.devnull, stderr=self.devnull)
            if not os.path.isfile('{}.hash'.format(linkfile)):
                call(sample.commands.fahash, shell=True, stdout=self.devnull, stderr=self.devnull)
            if not os.path.isfile('{}.txt'.format(linkfile)):
                call(sample.commands.epcr, shell=True, stdout=self.devnull, stderr=self.devnull)
            self.epcrqueue.task_done()

    def epcrparse(self):
        """
        Parse the ePCR text file outputs
        """
        printtime('Parsing ePCR results', self.start)
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'stx' in sample.general.datastore:
                    # Initialise count - this allows for the population of vtyperresults with unique values
                    uniquecount = 0
                    # This populates vtyperresults with the verotoxin subtypes
                    toxinlist = []
                    if os.path.isfile(sample[self.analysistype].resultsfile):
                        epcrresults = open(sample[self.analysistype].resultsfile, 'r')
                        for result in epcrresults:
                            # Only the lines without a # contain results
                            if "#" not in result:
                                uniquecount += 1
                                # Split on \t
                                data = result.split('\t')
                                # The subtyping primer pair is the first entry on lines with results
                                vttype = data[0].split('_')[0]
                                # Push the name of the primer pair - stripped of anything after a _ to the dictionary
                                if vttype not in toxinlist:
                                    toxinlist.append(vttype)

                    # Create a string of the entries in list1 joined with ";"
                    toxinstring = ";".join(sorted(toxinlist))
                    # Save the string to the metadata
                    sample[self.analysistype].toxinprofile = toxinstring
                else:
                    setattr(sample, self.analysistype, GenObject())
                    sample[self.analysistype].toxinprofile = 'NA'
            else:
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].toxinprofile = 'NA'

    def __init__(self, inputobject, analysistype):
        from Queue import Queue
        import multiprocessing
        self.metadata = inputobject.runmetadata.samples
        self.analysistype = analysistype
        self.reffilepath = inputobject.reffilepath
        self.start = inputobject.starttime
        # If the reference file path has not been provided, use the primer file supplied in the calling script
        # This is to allow the pipeline to run with a default primer file, and the stand-alone version to use
        # a custom primer file
        if not self.reffilepath:
            self.primerfile = inputobject.primerfile
        self.cpus = int(multiprocessing.cpu_count())
        self.epcrqueue = Queue(maxsize=self.cpus)
        self.devnull = open(os.devnull, 'wb')
        self.vtyper()


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':

    class Setup(object):

        def setup(self):
            """
            Set up the metadata object to be passed to Vtyper()
            """
            from glob import glob
            files = sorted(glob('{}*.fasta'.format(self.sequencepath)))
            samples = list()
            # Create the metadata for each file
            for fasta in files:
                # Create a metadata object to store all metadata associated with each strain
                metadata = MetadataObject()
                metadata.general = GenObject()
                metadata.commands = GenObject()
                # Set the name
                metadata.name = os.path.basename(fasta).split('.')[0]
                metadata.general.bestassemblyfile = fasta
                metadata.general.stx = True
                metadata.general.outputdirectory = self.path
                metadata.general.filenoext = fasta.split('.')[0]
                metadata.general.fastqfiles = list()
                samples.append(metadata)
            return samples

        def reporter(self):
            """
            Create a report of the results
            """
            printtime('Writing report', self.starttime)
            data = 'Strain,Profile\n'
            for sample in self.runmetadata.samples:
                # Only add to the string if there are results
                if sample[self.analysistype].toxinprofile:
                    data += '{},{}\n'.format(sample.name, sample[self.analysistype].toxinprofile)
            # Create the report, and write to it
            with open('{}/{}.csv'.format(self.reportpath, self.analysistype), 'wb') as report:
                report.write(data)

        def __init__(self):
            from argparse import ArgumentParser
            from time import time
            # Parser for arguments
            parser = ArgumentParser(
                description='Performs ePCR using a supplied primer file. The primers must be in the format: '
                            '<name>\t<forward primer>\t<reverse primer>\t<max size allowed between primers>\n.'
                            'Sequence files must be stored in <path>/sequences'
            )
            parser.add_argument('path',
                                help='Specify path in which reports are to be stored')
            parser.add_argument('-s', '--sequencepath',
                                required=True,
                                help='Path to assembly files')
            parser.add_argument('-f', '--primerfile',
                                required=True,
                                help='The name and path of the file containing the primers')
            # Get the arguments into an object
            arguments = parser.parse_args()
            self.starttime = time()
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(arguments.path, '')
            self.sequencepath = os.path.join(arguments.sequencepath, '')
            self.primerfile = arguments.primerfile
            # Initialise variables
            self.runmetadata = MetadataObject()
            self.reffilepath = False
            self.analysistype = 'ePCR'
            self.reportpath = os.path.join(self.path, 'reports')
            make_path(self.reportpath)
            # Initialise metadata
            self.runmetadata.samples = self.setup()
            # Run the analyses
            Vtyper(self, self.analysistype)
            # Create a report
            self.reporter()
            # Print the metadata to file
            printtime('Printing metadata to file', self.starttime)
            metadataprinter.MetadataPrinter(self)
            # Print a bold, green exit statement
            print u'\033[92m' + u'\033[1m' + u'\nElapsed Time: %0.2f seconds' % (time() - self.starttime) + u'\033[0m'

    Setup()

# -f /home/blais/git/V-Typer/vtx_subtyping_primers.txt
