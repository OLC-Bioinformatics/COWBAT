#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import printtime, make_path
from glob import glob
import shutil
from Bio import SeqIO
import time
import os
__author__ = 'adamkoziol'


class Get(object):

    def getrmlsthelper(self):
        """
        Makes a system call to rest_auth.pl, a Perl script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """
        from subprocess import call
        printtime('Downloading {} alleles'.format(self.analysistype), self.start)
        # Extract the path of the current script from the full path + file name
        homepath = os.path.split(os.path.abspath(__file__))[0]
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = os.path.join(self.path, self.analysistype)
        # System call
        rmlstupdatecall = 'cd {} && perl {} -a {}'\
            .format(newfolder, os.path.join(homepath, 'rest_auth.pl'), os.path.join(homepath, 'secret.txt'))
        # Create the path
        make_path(newfolder)
        # Copy over the access token to be used in the authentication
        shutil.copyfile(os.path.join(homepath, 'access_token'), os.path.join(newfolder, 'access_token'))
        # Run rest_auth.pl
        call(rmlstupdatecall, shell=True)
        # Get the new alleles into a list, and create the combinedAlleles file
        alleles = glob(os.path.join(newfolder, '*.tfa'))
        self.combinealleles(newfolder, alleles)

    def combinealleles(self, allelepath, alleles):
        printtime('Creating combined rMLST allele file', self.start)
        with open('{}/rMLST_combined.fasta'.format(allelepath), 'w') as combinedfile:
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
    from argparse import ArgumentParser
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
