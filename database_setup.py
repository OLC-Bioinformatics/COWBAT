#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import clear_logfile, combinetargets, MetadataObject, make_path, printtime, \
    run_subprocess, write_to_logfile
import get.get_rmlst as get_rmlst
import get.get_mlst as get_mlst
from argparse import ArgumentParser
from time import time
from glob import glob
import fileinput
import tarfile
import shutil
import os
__author__ = 'adamkoziol'


class DatabaseSetup(object):

    def main(self):
        """
        Run the methods
        """
        self.olc_databases()
        self.confindr()
        self.plasmidextractor()
        self.clark()
        self.mash()
        self.rmlst()
        self.mlst()
        self.cge_db_downloader('plasmidfinder', 'plasmidfinder_db', 'fsa', 'tfa')
        self.cge_db_downloader('resfinder', 'resfinder_db', 'fsa', 'tfa')
        # self.notes()
        self.cge_db_downloader('virulence', 'virulencefinder_db', 'fsa', 'tfa')
        self.cge_db_downloader('serosippr', 'serotypefinder_db', 'fsa', 'tfa')
        self.univec()

    def olc_databases(self):
        """
        Clone the OLC-specific databases from github. This method must be performed first, as the call will only clone
        the repository into an empty folder
        """
        printtime('Downloading OLC databases', self.start)
        # Set the git clone system call
        targetcall = 'git clone https://github.com/OLC-Bioinformatics/Databases.git {dbpath}'\
            .format(dbpath=self.databasepath)
        # Download the databases
        self.database_download(targetcall, self.databasepath)
        # Extract the databases from the archives
        printtime('Extracting databases from archives', self.start)
        for gz in glob(os.path.join(self.databasepath, '*.gz')):
            with tarfile.open(gz, 'r:gz') as tar:
                # Decompress the archive
                tar.extractall(path=self.databasepath)
            # Delete the archive file
            os.remove(gz)

    def confindr(self):
        """
        Download and extract the .tar file containing the ConFindr database from figshare.
        """
        databasepath = self.create_database_folder('ConFindr')
        tar_file = os.path.join(databasepath, 'confindr.tar')
        targetcall = 'wget -O {out} https://ndownloader.figshare.com/files/11864267'\
            .format(out=tar_file)
        self.database_download(targetcall, databasepath, complete=False)
        # Extract the databases from the archives
        printtime('Extracting database from archives', self.start)
        if not os.path.isfile(os.path.join(databasepath, 'complete')):
            with tarfile.open(tar_file, 'r') as tar:
                # Decompress the archive
                tar.extractall(path=databasepath)
            # Delete the archive file
            os.remove(tar_file)

    def plasmidextractor(self):
        """
        Download and extract the .tar file containing the PlasmidExtractor database from figshare.
        """
        databasepath = self.create_database_folder('plasmidextractor')
        tar_file = os.path.join(databasepath, 'confindr.tar')
        targetcall = 'wget -O {out} https://ndownloader.figshare.com/files/9827323'\
            .format(out=tar_file)
        self.database_download(targetcall, databasepath, complete=False)
        # Extract the databases from the archives
        printtime('Extracting database from archives', self.start)
        if not os.path.isfile(os.path.join(databasepath, 'complete')):
            with tarfile.open(tar_file, 'r') as tar:
                # Decompress the archive
                tar.extractall(path=databasepath)
            # Delete the archive file
            os.remove(tar_file)

    def clark(self):
        """
        Download and set-up the CLARK database using the set_targets.sh script. Use defaults of bacteria for database
        type, and species for taxonomic level
        """
        # Create the folder in which the database is to be stored
        databasepath = self.create_database_folder('clark')
        # Set the call to create the database - use the --light option, as we don't require the full database
        targetcall = 'cd {clarkpath} && ../opt/clark/set_targets.sh {dbpath} bacteria --species --light'\
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
        # Download the assembly summary refseq document
        summarycall = 'curl -o {} ftp://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'\
            .format(os.path.join(databasepath, 'assembly_summary_refseq.txt'))
        self.database_download(summarycall, databasepath, False)
        # Set the call to create the database
        targetcall = 'curl -o {} https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh'\
            .format(os.path.join(databasepath, 'RefSeqSketchesDefaults.msh'))
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

    def mlst(self, genera={'Escherichia', 'Vibrio', 'Campylobacter', 'Listeria', 'Bacillus', 'Staphylococcus',
                           'Salmonella'}):
        """
        Download the necessary up-to-date MLST profiles and alleles
        """
        printtime('Downloading MLST databases', self.start)
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

    def cge_db_downloader(self, analysistype, dbname, extension_in, extension_out):
        """
        Clones CGE databases into appropriate folder. Creates properly formatted file with non-redundant sequences
        :param analysistype: The name of the database folder to create
        :param dbname: The name of the database repository on bitbucket
        :param extension_in: The file extension of the FASTA files in the database
        :param extension_out: The desired extension for the FASTA files
        """
        printtime('Downloading {} database'.format(analysistype), self.start)
        if analysistype == 'serosippr':
            databasepath = os.path.join(self.databasepath, analysistype, 'Escherichia')
        else:
            databasepath = os.path.join(self.databasepath, analysistype)
        targetcall = 'git clone https://bitbucket.org/genomicepidemiology/{db}.git {atype}'\
            .format(db=dbname,
                    atype=databasepath)
        # Download the database
        self.database_download(targetcall, databasepath)
        # Create a variable to use in creating the combined targets file
        extension = extension_in
        # If the extension_out is different than extension_in, rename the files to have the appropriate extension
        if extension_in != extension_out:
            # Create a list of all the FASTA files with the input extension
            fastafiles = glob(os.path.join(databasepath, '*.{ext}'.format(ext=extension_in)))
            for fasta in fastafiles:
                # Split the extension
                filename = os.path.splitext(fasta)[0]
                # Rename the files
                os.rename(fasta, '{fn}.{ex}'.format(fn=filename,
                                                    ex=extension_out))
            # Update the variable to use when creating the combined targets file
            extension = extension_out
        # Create the combined targets file to use in the OLC pipelines
        if not os.path.isfile(os.path.join(databasepath, 'combinedtargets.fasta')):
            # Create the combinedtargets.fasta file - this will combine all the FASTA files in the downloaded database
            # into a properly-formatted, non-redundant FASTA database
            databasefiles = glob(os.path.join(databasepath, '*.{ext}'.format(ext=extension)))
            combinetargets(databasefiles, databasepath)

    def notes(self):
        """
        Clean the notes.txt file that comes with the resfinder database; it contains certain definitions with commas.
        Replace all commas in the file with semicolons. Append necessary changes to end of file
        """
        databasepath = self.create_database_folder('resfinder')
        # Load the file containing necessary changes to the resfinder notes.txt file
        with open(os.path.join(self.databasepath, 'resfinder_changes.txt'), 'r') as res_changes:
            changes = res_changes.read().splitlines()
        # Create a list to store all the entries in the notes file - will be used to ensure that the changes are not
        # entered more than once
        lines = list()
        # Open the notes file with fileinput to allow inplace editing
        note_file = os.path.join(databasepath, 'notes.txt')
        with fileinput.FileInput(note_file, inplace=True, backup='.bak') \
                as notes:
            for line in notes:
                # Replace all commas with semicolons
                print(line.replace(',', ';'), end='')
                lines.append(line)
        # Append changes to the end of the notes file
        with open(note_file, 'a') as notes:
            for change in changes:
                # Ensure that the changes have not already been added once
                if change + '\n' not in lines:
                    notes.write('{change}\n'.format(change=change))

    def univec(self):
        """
        Download the UniVec core database using wget
        """
        databasepath = self.create_database_folder('univec')
        #
        outputfile = os.path.join(databasepath, 'UniVec_core.tfa')
        targetcall = 'wget -O {} ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'\
            .format(outputfile)
        self.database_download(targetcall, databasepath)
        # Create a copy of the file with a .fasta extension
        if os.path.isfile(outputfile):
            renamed = os.path.splitext(outputfile)[0] + '.fasta'
            shutil.copy(outputfile, renamed)

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

    def database_download(self, targetcall, databasepath, complete=True):
        """
        Checks to see if the database has already been downloaded. If not, downloads the database, and writes stdout
        and stderr to the logfile
        :param targetcall: system call to download, and possibly set-up the database
        :param databasepath: absolute path of the database
        :param complete: boolean variable to determine whether the complete file should be created
        """
        # Create a file to store the logs; it will be used to determine if the database was downloaded and set-up
        completefile = os.path.join(databasepath, 'complete')
        # Run the system call if the database is not already downloaded
        if not os.path.isfile(completefile):
            out, err = run_subprocess(targetcall)
            print(out, err)
            # Write the out and err streams to the master files
            write_to_logfile(out, err, self.logfile, None, None, None, None)
            if complete:
                # Create the database completeness assessment file and populate it with the out and err streams
                with open(completefile, 'w') as complete:
                    complete.write(out)
                    complete.write(err)

    def __init__(self, args):
        self.databasepath = os.path.join(args.databasepath)
        make_path(self.databasepath)
        self.start = args.start
        # Determine the location of the CLARK scripts
        self.clarkpath = os.path.dirname(shutil.which('CLARK'))
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