#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import printtime, make_path
from Bio import SeqIO
from argparse import ArgumentParser
from glob import glob
import time
import os

from . import rest_auth_class
__author__ = 'adamkoziol'


class Get(object):

    def getrmlsthelper(self):
        """
        Makes a system call to rest_auth.py, a Python script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """

        printtime('Downloading {} alleles'.format(self.analysistype), self.start)
        # Extract the path of the current script from the full path + file name
        homepath = os.path.split(os.path.abspath(__file__))[0]
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = os.path.join(self.path, self.analysistype)
        # Create the path
        make_path(newfolder)
        # Create arguments to feed into the rest_auth_class script
        args = ArgumentParser
        args.secret_file = os.path.join(homepath, 'secret.txt')
        args.file_path = homepath
        args.output_path = newfolder
        args.start = self.start
        rmlst = rest_auth_class.REST(args)
        # Download the profile and alleles
        rmlst.main()

        # Get the new alleles into a list, and create the combinedAlleles file
        alleles = glob(os.path.join(newfolder, '*.tfa'))
        self.combinealleles(newfolder, alleles)

    def combinealleles(self, allelepath, alleles):
        printtime('Creating combined rMLST allele file', self.start)
        with open(os.path.join(allelepath, 'rMLST_combined.fasta'), 'w') as combinedfile:
            # Open each allele file
            for allele in sorted(alleles):
                # with open(allele, 'rU') as fasta:
                for record in SeqIO.parse(open(allele, "rU"), "fasta"):
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
                    # Write each record to the combined file
                    SeqIO.write(record, combinedfile, 'fasta')

    def __init__(self, args):
        self.path = os.path.join(args.path)
        self.start = args.start
        self.analysistype = 'rMLST'
        self.getrmlsthelper()


if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    # Parser for arguments
    parser = ArgumentParser(description='')
    parser.add_argument('path',
                        help='Specify input directory')

    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()

    # Run the script
    Get(arguments)

    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')
