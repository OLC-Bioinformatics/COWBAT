#!/usr/bin/env python
from subprocess import Popen, PIPE, STDOUT
__author__ = 'adamkoziol'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    import os
    import errno
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    from collections import defaultdict
    return defaultdict(make_dict)


def printtime(string, start):
    """Prints a string in bold with the elapsed time
    :param string: a string to be printed in bold
    :param start: integer of the starting time
    """
    import time
    print('\n\033[1m' + "[Elapsed Time: {:.2f} seconds] {}".format(time.time() - start, string) + '\033[0m')

# Initialise globalcount
globalcount = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    import time
    import sys
    # Use a global variable
    global globalcount
    if globalcount <= 80:
        sys.stdout.write('.')
        globalcount += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        globalcount = 1


def execute(command, outfile=""):
    """
    Allows for dots to be printed to the terminal while waiting for a long system call to run
    :param command: the command to be executed
    :param outfile: optional string of an output file
    from https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
    """
    import sys
    import time
    # Initialise count
    count = 0
    # Initialise the starting time
    start = int(time.time())
    maxtime = 0
    print command
    # Run the commands - direct stdout to PIPE and stderr to stdout
    process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
    # Write the initial time
    sys.stdout.write('[{:}] '.format(time.strftime('%H:%M:%S')))
    # Create the output file - if not provided, then nothing should happen
    writeout = open(outfile, "ab+") if outfile else ""
    # Poll process for new output until finished
    while True:
        # If an output file name is provided
        if outfile:
            # Get stdout into a variable
            nextline = process.stdout.readline()
            # Print stdout to the file
            writeout.write(nextline)
        # Break from the loop if the command is finished
        if process.poll() is not None:
            break
        # Adding sleep commands slowed down this method when there was lots of output. Difference between the start time
        # of the analysis and the current time. Action on each second passed
        currenttime = int(time.time())
        if currenttime - start > maxtime:
            # Set the max time for each iteration
            maxtime = currenttime - start
            # Print up to 80 dots on a line, with a one second delay between each dot
            if count <= 80:
                sys.stdout.write('.')
                count += 1
            # Once there are 80 dots on a line, start a new line with the the time
            else:
                sys.stdout.write('\n[{:}] .'.format(time.strftime('%H:%M:%S')))
                count = 1
    # Close the output file
    writeout.close() if outfile else ""
    sys.stdout.write('\n')


def filer(filelist):
    """
    Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)
    :param filelist: List of files to parse
    """
    import re
    # Initialise the set
    fileset = set()
    for seqfile in filelist:
        # Search for the conventional motifs present following strain names
        # _S\d+_L001_R\d_001.fastq(.gz) is a typical unprocessed Illumina fastq file
        if re.search("_S\d+_L001", seqfile):
            fileset.add(re.split("_S\d+_L001", seqfile)[0])
        # Files with _R\d_001.fastq(.gz) are created in the SPAdes assembly pipeline
        elif re.search("_R\d_001", seqfile):
            fileset.add(re.split("_R\d_001", seqfile)[0])
        # _R\d.fastq(.gz) represents a simple naming scheme for paired end reads
        elif re.search("R\d.fastq", seqfile):
            fileset.add(re.split("_R\d.fastq", seqfile)[0])
        # _\d.fastq is always possible
        elif re.search("[-_]\d.fastq", seqfile):
            fileset.add(re.split("[-_]\d.fastq", seqfile)[0])
        # .fastq is the last option
        else:
            fileset.add(re.split(".fastq", seqfile)[0])
        dotter()
    return fileset


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = (lambda y: y if y else {})(x)
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        if value:
            self.datastore[key] = value
        else:
            self.datastore[key] = "NA"


class MetadataObject(object):
    """Object to store static variables"""
    def __init__(self):
        """Create datastore attr with empty dict"""
        super(MetadataObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        """:key is retrieved from datastore if exists, for nested attr recursively :self.__setattr__"""
        if key not in self.datastore:
            self.__setattr__(key)
        return self.datastore[key]

    def __setattr__(self, key, value=GenObject(), **args):
        """Add :value to :key in datastore or create GenObject for nested attr"""
        if args:
            self.datastore[key].value = args
        else:
            self.datastore[key] = value

    def dump(self):
        """Print only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = {}
        for attr in self.datastore:
            metadata[attr] = {}
            if not attr.startswith('__'):
                if isinstance(self.datastore[attr], str):
                    metadata[attr] = self.datastore[attr]
                else:
                    metadata[attr] = self.datastore[attr].datastore
        return metadata
