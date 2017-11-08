#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import clear_logfile, combinetargets, MetadataObject, make_path, printtime, \
    run_subprocess, write_to_logfile
import get.get_rmlst as get_rmlst
import get.get_mlst as get_mlst
from argparse import ArgumentParser
from time import time
from glob import glob
import shutil
import os
__author__ = 'adamkoziol'


class DatabaseSetup(object):

    def main(self):
        """
        Run the methods
        """
        self.clark()
        self.mash()
        self.rmlst()
        self.mlst()
        self.cge_db_downloader('plasmidfinder', 'plasmidfinder_db', 'fsa')
        self.cge_db_downloader('resfinder', 'resfinder_db', 'fsa')
        self.cge_db_downloader('virulence', 'virulencefinder_db', 'fsa')
        self.cge_db_downloader('serosippr', 'serotypefinder_db', 'fsa')
        self.univec()

    def clark(self):
        """
        Download and set-up the CLARK database using the set_targets.sh script. Use defaults of bacteria for database
        type, and species for taxonomic level
        """
        # Create the folder in which the database is to be stored
        databasepath = self.create_database_folder('clark')
        # Set the call to create the database
        targetcall = 'cd {clarkpath} && ./set_targets.sh {dbpath} bacteria --species'\
            .format(clarkpath=self.clarkpath,
                    dbpath=databasepath)
        # Download the database
        self.database_download(targetcall, databasepath)

    def mash(self):
        """
        Download the pre-computed sketch of the RefSeq database, and compress it with gzip
        """
        # Create the folder in which the database is to be stored
        databasepath = self.create_database_folder('mash')
        # Set the call to create the database
        targetcall = 'curl https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh | gzip > {}'\
            .format(os.path.join(databasepath, 'RefSeqSketchesDefault.msh.gz'))
        # Download the database
        self.database_download(targetcall, databasepath)

    def rmlst(self):
        """
        Get the most up-to-date profiles and alleles from pubmlst. Note that you will need the necessary access token
        and secret for this to work
        """
        printtime('Downloading rMLST database', self.start)
        # Set the name of the file to be used to determine if the database download and set-up was successful
        completefile = os.path.join(self.databasepath, 'rMLST', 'complete')
        if not os.path.isfile(completefile):
            # Create an object to send to the rMLST download script
            args = MetadataObject()
            # Add the path and start time attributes
            args.path = self.databasepath
            args.start = self.start
            # Run the rMLST download
            get_rmlst.Get(args)
            # Create and populate the complete.txt file
            with open(completefile, 'w') as complete:
                complete.write('\n'.join(glob(os.path.join(self.databasepath, 'rMLST', '*'))))

    def mlst(self):
        """
        Download the necessary up-to-date MLST profiles and alleles
        """
        printtime('Downloading MLST databases', self.start)
        genera = ['Escherichia', 'Vibrio', 'Campylobacter', 'Listeria', 'Bacillus']
        for genus in genera:
            # Create an object to pass to the get_mlst script
            args = MetadataObject()
            # Populate the object with the necessary attributes
            args.species = genus
            args.repository_url = 'http://pubmlst.org/data/dbases.xml'
            args.force_scheme_name = False
            args.path = os.path.join(self.databasepath, 'MLST', genus)
            # Create the name of the file to be used to determine if the database download and setup was successful
            completefile = os.path.join(args.path, 'complete')
            # Only download the files if the download was not previously successful
            if not os.path.isfile(completefile):
                # Run the download
                get_mlst.main(args)
                # Create and populate the complete.txt file
                with open(completefile, 'w') as complete:
                    complete.write('\n'.join(glob(os.path.join(args.path, '*'))))

    def cge_db_downloader(self, analysistype, dbname, extension):
        """
        Clones CGE databases into appropriate folder. Creates properly formatted file with non-redundant sequences
        :param analysistype: The name of the database folder to create
        :param dbname: The name of the database repository on bitbucket
        :param extension: The file extension of the FASTA files
        """
        printtime('Downloading {} database'.format(analysistype), self.start)
        databasepath = os.path.join(self.databasepath, analysistype)
        targetcall = 'git clone https://bitbucket.org/genomicepidemiology/{db}.git {atype}'\
            .format(db=dbname,
                    atype=databasepath)
        # Download the database
        self.database_download(targetcall, databasepath)
        if not os.path.isfile(os.path.join(databasepath, 'combinedtargets.fasta')):
            # Create the combinedtargets.fasta file - this will combine all the FASTA files in the downloaded database
            # into a properly-formatted, non-redundant FASTA database
            databasefiles = glob(os.path.join(databasepath, '*.{ext}'.format(ext=extension)))
            combinetargets(databasefiles, databasepath)

    def univec(self):
        """
        Download the UniVec core database using wget
        """
        databasepath = self.create_database_folder('univec')
        targetcall = 'wget -O {} ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'\
            .format(os.path.join(databasepath, 'UniVec_core.tfa'))
        self.database_download(targetcall, databasepath)

    def create_database_folder(self, database):
        """
        Create an appropriately named folder in which the database is to be stored
        :param database: the name of the database folder to create
        :return: the absolute path of the folder
        """
        printtime('Setting up {} database'.format(database), self.start)
        # Define the path to store the database files
        databasepath = os.path.join(self.databasepath, database)
        # Create the path as required
        make_path(databasepath)
        return databasepath

    def database_download(self, targetcall, databasepath):
        """
        Checks to see if the database has already been downloaded. If not, downloads the database, and writes stdout
        and stderr to the logfile
        :param targetcall: system call to download, and possibly set-up the database
        :param databasepath: absolute path of the database
        """
        # Create a file to store the logs; it will be used to determine if the database was downloaded and set-up
        completefile = os.path.join(databasepath, 'complete')
        # Run the system call if the database is not already downloaded
        if not os.path.isfile(completefile):
            out, err = run_subprocess(targetcall)
            # Write the out and err streams to the master files
            write_to_logfile(out, err, self.logfile, None, None, None, None)
            # Create the database completeness assessment file and populate it with the out and err streams
            with open(completefile, 'w') as complete:
                complete.write(out)
                complete.write(err)

    def __init__(self, args):
        self.databasepath = os.path.join(args.databasepath)
        make_path(self.databasepath)
        self.start = args.start
        # Determine the location of the CLARK scripts
        self.clarkpath = os.path.dirname(shutil.which('estimate_abundance.sh'))
        self.logfile = os.path.join(self.databasepath, 'logfile')
        # Delete log files form previous iterations of the script in this folder
        clear_logfile(self.logfile)


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Downloads and sets up required databases for the CFIA OLC Workflow '
                                        'for Bacterial Assembly and Typing (COWBAT)')
    parser.add_argument('-d', '--databasepath',
                        required=True,
                        help='Absolute path to location to store database files. Include any version numbers if '
                             'required.')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.start = time()
    # Run the pipeline
    pipeline = DatabaseSetup(arguments)
    pipeline.main()
    printtime('Databases downloaded and set up', arguments.start)
