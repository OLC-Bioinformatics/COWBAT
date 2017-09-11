#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import *
from glob import glob
import shutil
from Bio import SeqIO
import time
__author__ = 'adamkoziol'


class Get(object):

    def getrmlsthelper(self, referencefilepath, update, start):
        """
        Makes a system call to rest_auth.pl, a Perl script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """
        from subprocess import call
        analysistype = 'rMLST'
        # Folders are named based on the download date e.g 2016-04-26
        # Find all folders (with the trailing / in the glob search) and remove the trailing /
        lastfolder = sorted(glob('{}{}/2*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
        delta, foldersize, d1 = self.schemedate(lastfolder)
        # Extract the path of the current script from the full path + file name
        homepath = os.path.split(os.path.abspath(__file__))[0]
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = '{}{}/{}'.format(referencefilepath, analysistype, d1)
        # System call
        rmlstupdatecall = 'cd {} && perl {}/rest_auth.pl -a {}/secret.txt'.format(newfolder, homepath, homepath)
        if update:
            if delta.days > 7 or foldersize < 100:
                printtime("Last update of rMLST profile and alleles was {} days ago. Updating".format(str(delta.days)),
                          start)
                # Create the path
                make_path(newfolder)
                # Copy over the access token to be used in the authentication
                shutil.copyfile('{}/access_token'.format(homepath), '{}/access_token'.format(newfolder))
                # Run rest_auth.pl
                call(rmlstupdatecall, shell=True)
                # Get the new alleles into a list, and create the combinedAlleles file
                alleles = glob('{}/*.tfa'.format(newfolder))
                self.combinealleles(start, newfolder, alleles)
            # If the profile and alleles are up-to-date, set :newfolder to :lastfolder
            else:
                newfolder = lastfolder
            # Ensure that the profile/alleles updated successfully
            # Calculate the size of the folder by adding the sizes of all the files within the folder together
            newfoldersize = sum(os.path.getsize('{}/{}'.format(newfolder, f)) for f in os.listdir(newfolder)
                                if os.path.isfile('{}/{}'.format(newfolder, f)))
            # If the profile/allele failed, remove the folder, and use the most recent update
            if newfoldersize < 100:
                shutil.rmtree(newfolder)
                newfolder = sorted(glob('{}{}/*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
        # Don't update the profile/alleles if not requested
        else:
            newfolder = lastfolder
        # Return the system call and the folder containing the profile and alleles
        return rmlstupdatecall, newfolder

    @staticmethod
    def schemedate(lastfolder):
        from datetime import date
        try:
            # Extract the folder name (date) from the path/name
            lastupdate = os.path.split(lastfolder)[1]
        except AttributeError:
            lastupdate = '2000-01-01'
        try:
            # Calculate the size of the folder by adding the sizes of all the files within the folder together
            foldersize = sum(os.path.getsize('{}/{}'.format(lastfolder, f)) for f in os.listdir(lastfolder)
                             if os.path.isfile('{}/{}'.format(lastfolder, f)))
        except TypeError:
            foldersize = 0
        # Try to figure out the year, month, and day from the folder name
        try:
            (year, month, day) = lastupdate.split("-")
            # Create a date object variable with the year, month, and day
            d0 = date(int(year), int(month), int(day))
        except ValueError:
            # Set an arbitrary date in the past to force an update
            d0 = date(2000, 1, 1)
        # Create a date object with the current date
        d1 = date(int(time.strftime("%Y")), int(time.strftime("%m")), int(time.strftime("%d")))
        # Subtract the last update date from the current date
        delta = d1 - d0

        return delta, foldersize, d1

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
