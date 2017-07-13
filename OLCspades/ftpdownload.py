#!/usr/bin/env python
from accessoryFunctions import *
from cStringIO import StringIO
import re
from time import sleep
import pycurl
__author__ = 'adamkoziol'


class Download(object):

    def main(self):
        """
        Run the necessary methods in the proper order
        """
        printtime('Finding sequencing runs on the FTP', self.starttime)
        self.findfiles()
        printtime('Downloading file(s) from FTP', self.starttime)
        self.ftp()
        printtime('Decompressing sequencing run(s)', self.starttime)
        self.decompress()
        printtime('Renaming folder(s) to allow automatic assembly', self.starttime)
        self.renamefolder()

    def findfiles(self):
        """
        Search for sequencing runs on the FTP
        """
        # Load the credentials to the FTP
        self.username, self.password = open(self.ftpcredentials, 'r').readline().rstrip().split(',')
        # Create a URL that includes the user name and password, so PycURL can login to the FTP server
        destinationurl = 'ftp://{}:{}@{}'.format(self.username, self.password, self.downloadpath)
        success = False
        while not success:
            # Create a pycurl object
            c = pycurl.Curl()
            # Specify the details of FTP server
            c.setopt(pycurl.URL, destinationurl)
            # Create a buffer to store the output
            output = StringIO()
            # Assign this buffer to pycurl object
            c.setopt(pycurl.WRITEFUNCTION, output.write)
            # Perform the LIST operation
            c.perform()
            c.close()
            # Get the output in a string
            result = output.getvalue()
            # FTP LIST output is separated by \n - split the output on newlines
            lines = result.split('\n')
            for line in lines:
                if line:
                    # Split the file name from the end of the list
                    zipfile = line.split()[-1]
                    if zipfile.endswith('.zip'):
                        # Only add .zip files to the list of files to download
                        self.files.append(os.path.join(self.downloadpath, zipfile))
                        success = True
            # Wait for an hour before checking again
            if not self.files:
                printtime('No runs currently on the FTP. Will check again in 60 minutes', self.starttime)
                sleep(3600)

    def ftp(self):
        """
        Uses PycURL to download the zip file containing the sequencing run from the FTP 
        """
        for zipfile in self.files:
            # Create an object to store metadata for each file
            metadata = MetadataObject()
            # Set the attributes for the object
            metadata.localfile = os.path.join(self.assemblypath, os.path.basename(zipfile))
            metadata.decompressed = metadata.localfile.split('.')[0]
            metadata.readyfolder = metadata.decompressed + '_Ready'
            metadata.extractpath = os.path.join(self.assemblypath, os.path.basename(zipfile).split('.')[0])
            # Count downloaded size
            count = 0
            downloading = True
            success = False
            while not success and downloading:
                try:
                    # Create a StringIO instance to store data from the curl request
                    header = StringIO()
                    # Create a pycurl object to determine the size of the file to download
                    curlcheck = pycurl.Curl()
                    # Use the .setopt attribute to set options
                    # The ftp link
                    curlcheck.setopt(curlcheck.URL, zipfile)
                    # Only get the headers - not the body
                    curlcheck.setopt(curlcheck.NOBODY, 1)
                    # Write the response to the StringIO instance
                    curlcheck.setopt(curlcheck.WRITEFUNCTION, header.write)
                    # Run the curl command
                    curlcheck.perform()
                    # Close
                    curlcheck.close()
                    try:
                        # Pull the filesize from the header information
                        filesize = int(re.findall('Content-Length: (.+)\r\n', header.getvalue())[0])
                        # If either the compressed or decompressed file is present, determine the filesize
                        if os.path.isfile(metadata.localfile) or os.path.isfile(metadata.decompressed):
                            count = os.path.getsize(metadata.localfile) if os.path.isfile(metadata.localfile) else \
                                os.path.getsize(metadata.decompressed)
                        # If the local file size is the same as the filesize on the ftp, then don't download
                        if count >= filesize:
                            # Already downloaded
                            downloading = False
                            success = True
                        else:
                            # Open the destination file to write
                            with open(metadata.localfile, 'wb') as localfile:
                                # Create a pycurl instance to download the file
                                curlinstance = pycurl.Curl()
                                # # Set the desired encoding type to be gzip
                                curlinstance.setopt(curlinstance.URL, zipfile)
                                # Write the data to the download file
                                curlinstance.setopt(curlinstance.WRITEDATA, localfile)
                                curlinstance.perform()
                                curlinstance.close()
                                # Once finished, set success to True to break the while loop
                                success = True
                                downloading = False
                    except IndexError as e:
                        print('\nFile does not exist, retrying in 10 minutes: ' + str(e))
                        sleep(600)
                except IOError as e:
                    print('\nDownload error, retrying in one minute: ' + str(e))
                    sleep(60)
            self.runmetadata.append(metadata)

    def decompress(self):
        """
        Decompress the contents of the zip file to a folder in the assembly path
        """
        import zipfile
        for sample in self.runmetadata:
            if not os.path.isdir(sample.readyfolder):
                # Open the compressed file, and decompress the contents to the extraction destination
                with zipfile.ZipFile(sample.localfile, "r") as zip_ref:
                    zip_ref.extractall(sample.extractpath)

    def renamefolder(self):
        """
        Rename the folder by appending _Ready to let the automated assembly pipeline wrapper
        know that this run is ready to be assembled
        """
        for sample in self.runmetadata:
            # Don't try to rename the folder if it has previously been renamed
            if not os.path.isdir(sample.readyfolder) and not os.path.isdir(sample.readyfolder + '_Queued'):
                os.rename(sample.extractpath, sample.readyfolder)

    def __init__(self, args):
        """
        :param args: object of arguments passed to the script
        Initialises the variables required for this class
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.downloadpath = os.path.join(args.downloadpath, '')
        self.assemblypath = os.path.join(args.assemblypath, '')
        # Define the start time
        self.starttime = args.starttime
        self.ftpcredentials = os.path.join(args.homepath, 'credentials')
        self.username = str()
        self.password = str()
        self.files = list()
        self.runmetadata = list()
        # Assertions to ensure that the provided variables are valid
        assert os.path.isdir(self.assemblypath), u'Supplied assembly path location is not a valid directory {0!r:s}' \
            .format(self.assemblypath)
        # Perform the download
        self.main()

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    import os
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Download sequencing runs from the FTP, and prep them for automatic assembly')
    parser.add_argument('-d', '--downloadpath',
                        required=True,
                        help='Path to a folder containing zip files with sequencing runs to download anonymously from '
                             'the FTP (e.g.ftp.agr.gc.ca/incoming/cfia-ak/')
    parser.add_argument('-a', '--assemblypath',
                        required=True,
                        help='Path to place the uncompressed folder containing the sequencing run')
    # Get the arguments into an object
    arguments = parser.parse_args()
    arguments.starttime = time()
    # Extract the path of the current script from the full path + file name
    arguments.homepath = os.path.split(os.path.abspath(__file__))[0]
    # Run the pipeline
    Download(arguments)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time() - arguments.starttime) + '\033[0m')
