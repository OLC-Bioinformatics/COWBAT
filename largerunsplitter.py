#!/usr/bin/env python
import os
__author__ = 'adamkoziol'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    import errno
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


class SplitRun(object):

    def splitrun(self):
        from glob import glob
        # Get the fastq files into a list
        fastqfiles = sorted(glob('{}*fastq*'.format(self.path)))
        # Create necessary variables
        currentname = ''
        currentfastq = ''
        listoffastq = []
        # Iterate through the list of fastq files
        for i, fastqfile in enumerate(fastqfiles):
            # Split the fastq files into groups of 16 (using 32 because there are two fastq files per sample)
            if i % 32 == 0 or i == len(fastqfiles) - 1 and i != 0:
                # Set the name of the previous and current fastq files - to be used in naming the folders
                previousname = currentname
                currentname = os.path.split(fastqfile)[1].split('_')[0]
                # if previousname:
                foldername = '{}_{}'.format(previousname, currentfastq)
                print foldername
                # Create the folder at the destination
                make_path('{}{}'.format(self.runpath, foldername))
                # Clear the list of fastq files
                listoffastq = []
            # Set the name of the current fastq file - will be used in naming the folders
            currentfastq = os.path.split(fastqfile)[1].split('_')[0]
            # Append the fastq file to the list
            listoffastq.append(fastqfile)

    def __init__(self, args):
        from time import time
        starttime = time()
        # Initialise the variables
        self.path = os.path.join(args.path, '')
        self.runpath = os.path.join(args.destination, '')
        # Split the run
        self.splitrun()
        # Print an exit statement
        print "\nElapsed Time: %s seconds" % (time() - starttime)

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Runs kSNPv3automate.py on all the subfolders in a directory')
    parser.add_argument('path',  help='Specify path')
    parser.add_argument('-d', '--destination', default='/nas/akoziol/WGS_Spades/',
                        help='Optional folder path to store .fastq files created using the fastqCreation module. '
                             '/nas/akoziol/WGS_Spades/')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SplitRun(arguments)
