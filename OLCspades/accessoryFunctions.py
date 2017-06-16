#!/usr/bin/env python
import errno
import os
from subprocess import Popen, PIPE, STDOUT

# noinspection PyProtectedMember
from Bio.Application import _Option, AbstractCommandline, _Switch
__author__ = 'adamkoziol'


def get_version(exe):
    """
    :param exe: :type list required
    """
    assert isinstance(exe, list)
    return Popen(exe, stdout=PIPE, stderr=STDOUT).stdout.read()


def logstr(*args):
    yield "{}\n".__add__("-".__mul__(60)).__mul__(len(args)).format(*args)


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
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


class Dotter(object):

    def globalcounter(self):
        """Resets the globalcount to 0"""
        self.globalcount = 0

    def dotter(self):
        """Prints formatted time to stdout at the start of a line, as well as a "."
        whenever the length of the line is equal or lesser than 80 "." long"""
        import sys
        if self.globalcount <= 80:
            sys.stdout.write('.')
            self.globalcount += 1
        else:
            sys.stdout.write('\n.')
            self.globalcount = 1

    def __init__(self):
        self.globalcount = 0

# Initialise globalcount
globalcount = 0


def globalcounter():
    """Resets the globalcount to 0"""
    global globalcount
    globalcount = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # import time
    import sys
    # Use a global variable
    global globalcount
    if globalcount <= 80:
        sys.stdout.write('.')
        globalcount += 1
    else:
        sys.stdout.write('\n.')
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
    # Removing Shell=True to prevent excess memory use thus shlex split if needed
    if type(command) is not list:
        import shlex
        command = shlex.split(command)
    # Run the commands - direct stdout to PIPE and stderr to stdout
    # DO NOT USE subprocess.PIPE if not writing it!
    if outfile:
        process = Popen(command, stdout=PIPE, stderr=STDOUT)
    else:
        devnull = open(os.devnull, 'wb')
        process = Popen(command, stdout=devnull, stderr=STDOUT)
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


def filer(filelist, extension='fastq'):
    """
    Helper script that creates a set of the stain names created by stripping off parts of the filename.
    Hopefully handles different naming conventions (e.g. 2015-SEQ-001_S1_L001_R1_001.fastq(.gz),
    2015-SEQ-001_R1_001.fastq.gz, 2015-SEQ-001_R1.fastq.gz, 2015-SEQ-001_1.fastq.gz, and 2015-SEQ-001_1.fastq.gz
    all become 2015-SEQ-001)
    :param filelist: List of files to parse
    :param extension: the file extension to use. Default value is 'fastq
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
        elif re.search("R\d.{}".format(extension), seqfile):
            fileset.add(re.split("_R\d.{}".format(extension), seqfile)[0])
        # _\d.fastq is always possible
        elif re.search("[-_]\d.{}".format(extension), seqfile):
            fileset.add(re.split("[-_]\d.{}".format(extension), seqfile)[0])
        # .fastq is the last option
        else:
            fileset.add(re.split(".{}".format(extension), seqfile)[0])
    return fileset


def relativesymlink(src_file, dest_file):
    """
    https://stackoverflow.com/questions/9793631/creating-a-relative-symlink-in-python-without-using-os-chdir
    :param src_file: the file to be linked
    :param dest_file: the path and filename to which the file is to be linked
    """
    # Perform relative symlinking
    try:
        print(os.path.relpath(src_file), os.path.relpath(dest_file))
        os.symlink(
            # Find the relative path for the source file and the destination file
            os.path.relpath(src_file),
            os.path.relpath(dest_file)
        )
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = x if x else {}
        # start = (lambda y: y if y else {})(x)
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        if self.datastore[key] or self.datastore[key] == 0 or self.datastore[key] is False or all(self.datastore[key]):
            return self.datastore[key]
        else:
            self.datastore[key] = 'NA'
            return self.datastore[key]

    def __setattr__(self, key, value):
        if value:
            self.datastore[key] = value
        elif type(value) != int:
            if value is False or all(value):
                self.datastore[key] = value
        else:
            if value >= 0:
                self.datastore[key] = 0
            else:
                self.datastore[key] = "NA"

    def __delattr__(self, key):
        del self.datastore[key]


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

    def __getitem__(self, item):
        return self.datastore[item]

    def dump(self):
        """Prints only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = {}
        for attr in self.datastore:
            metadata[attr] = {}
            if not attr.startswith('__'):
                if isinstance(self.datastore[attr], str):
                    metadata[attr] = self.datastore[attr]
                else:
                    try:
                        metadata[attr] = self.datastore[attr].datastore
                    except AttributeError:
                        print(attr)
        return metadata


class MakeBlastDB(AbstractCommandline):
    """Base makeblastdb wrapper"""
    def __init__(self, cmd='makeblastdb', **kwargs):
        assert cmd is not None
        extra_parameters = [
            # Core:
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["-help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["-version", "version"],
                    "Print version number;  ignore other arguments."),
            # Output configuration options
            _Option(["-out", "out"],
                    "Output file prefix for db.",
                    filename=True,
                    equate=False),
            _Option(["-in", "db"],
                    "The sequence create db with.",
                    filename=True,
                    equate=False),  # Should this be required?
            _Option(["-dbtype", "dbtype"],
                    "Molecule type of target db (string, 'nucl' or 'prot').",
                    equate=False)]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


def combinetargets(targets, targetpath):
    """
    :param targets: fasta gene targets to combine
    :param targetpath: folder containing the targets
    """
    from Bio import SeqIO
    make_path(targetpath)
    recordlist = []
    # Open the combined allele file to write
    with open('{}/combinedtargets.tfa'.format(targetpath), 'w') as combinedfile:
        # Open each target file
        for target in sorted(targets):
            # with open(allele, 'rU') as fasta:
            for record in SeqIO.parse(open(target, 'rU'), 'fasta'):
                # Extract the sequence record from each entry in the multifasta
                # Replace and dashes in the record.id with underscores
                record.id = record.id.replace('-', '_')
                # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                # with gaps
                # noinspection PyProtectedMember
                record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                # Clear the name and description attributes of the record
                record.name = ''
                record.description = ''
                if record.id not in recordlist:
                    # Write each record to the combined file
                    SeqIO.write(record, combinedfile, 'fasta')
                    recordlist.append(record.id)


class KeyboardInterruptError(Exception):
    pass
