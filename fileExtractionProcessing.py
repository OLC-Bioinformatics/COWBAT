__author__ = 'akoziol'

# OS commands
import os
# Subprocess->call is used for making system calls
import subprocess
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob
# Shutil is useful for file moving/copying
import shutil
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# System tools
import sys


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def unzipping():
    zipFile = glob.glob("*.zip")
    if not zipFile:
        print("Files are already extracted from archive.")
        pass
    else:
        zipCommand = "unzip -j -qq %s" % zipFile[0]
        # print(zipCommand)
        print "Extracting files from archive"
        subprocess.call(zipCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # Remove the zipfile
        os.remove(zipFile[0])
        # Remove the "Undetermined" reads
        undeterminedFiles = glob.glob("Undetermined*")
        for file in undeterminedFiles:
            os.remove(file)


def moveExtract(strain, gzFiles, path):
    """Renames, moves, and uncompresses .gz files"""
    forward = str(strain) + "_R1_001.fastq.gz"
    reverse = str(strain) + "_R2_001.fastq.gz"
    shutil.move(gzFiles[0], "%s/%s/%s" % (path, strain, forward))
    shutil.move(gzFiles[1], "%s/%s/%s" % (path, strain, reverse))
    sys.stdout.write('.')
    gzipCommandForward = "gzip -d %s/%s/%s" % (path, strain, forward)
    gzipCommandReverse = "gzip -d %s/%s/%s" % (path, strain, reverse)
    subprocess.call(gzipCommandForward, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    subprocess.call(gzipCommandReverse, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))


def folderer(sampleNames, path):
    """Move the .gz files into appropriately named folders, and decompress them"""
    gzCheck = glob.glob("*.gz")
    if not gzCheck:
        print("Processing files.")
        for strain in sampleNames:
            os.chdir("%s/%s" % (path, strain))
            gzFiles = glob.glob("*.gz")
            if not gzFiles:
                pass
            else:
                moveExtract(strain, gzFiles, path)
    else:
        print("Moving and extracting fastq files.")
        for strain in sampleNames:
            # Make the required folders (if necessary)
            make_path(strain)
            # Get the .gz files into a list
            gzFiles = glob.glob("%s_S*" % strain)
            if not gzFiles:
                os.chdir("%s/%s" % (path, strain))
                gzFiles = glob.glob("*.gz")
                if not gzFiles:
                    pass
                else:
                    moveExtract(strain, gzFiles, path)
            else:
                moveExtract(strain, gzFiles, path)


def functionsGoNOW(sampleNames, path):
    """Run the functions"""
    sampleNames = sampleNames
    path = path
    unzipping()
    folderer(sampleNames, path)
