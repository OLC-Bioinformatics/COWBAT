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
from multiprocessing import Pool
import re

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def unzipping():
    """Unzips the zip archive of .fastq.gz files if it exists"""
    zipFile = glob.glob("*.zip")
    if not zipFile:
        print("Files are already extracted from archive.")
        undeterminedFiles = glob.glob("Undetermined*")
        for file in undeterminedFiles:
            os.remove(file)
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


def foldererPrepProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    foldererPrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'fileExtractionProcessing':
        createfoldererPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            foldererPrepArgs.append((name, path))
        # This map function allows for multi-processing
        createfoldererPool.map(folderer, foldererPrepArgs)


def moveExtract(strain, gzFiles, path, newPath, seqNum):
    """Renames, moves, and uncompresses .gz files"""
    forward = str(strain) + "_R1_001.fastq"
    reverse = str(strain) + "_R2_001.fastq"

    # uncompressedFiles = glob.glob("%s/%s/*.fastq" % (path, strain))
    if seqNum:
        # print strain, gzFiles[0]
    # if os.path.isfile("%s/%s" % (path, gzFiles[0])) and os.path.isfile("%s/%s" % (path, gzFiles[1])) and len(uncompressedFiles) < 2:
    # if gzFiles[0]:
        shutil.move("%s/%s" % (newPath, gzFiles[0]), "%s/%s/%s.gz" % (path, strain, forward))
    # if gzFiles[1]:
        shutil.move("%s/%s" % (newPath, gzFiles[1]), "%s/%s/%s.gz" % (path, strain, reverse))
    #     sys.stdout.write('.')
    if not os.path.isfile("%s/%s/%s" % (path, strain, forward)):
        gzipCommandForward = "gzip -d --force %s/%s/%s.gz" % (path, strain, forward)
        subprocess.call(gzipCommandForward, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    if not os.path.isfile("%s/%s/%s" % (path, strain, reverse)):
        gzipCommandReverse = "gzip -d --force %s/%s/%s.gz" % (path, strain, reverse)
        subprocess.call(gzipCommandReverse, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))


def folderer((name, path)):
    """Move the .gz files into appropriately named folders, and decompress them"""
    # os.chdir(path)
    # print name
    # os.system("pwd")
    gzCheck = [f for f in os.listdir(path) if re.search("%s\w+.fastq.gz" % name, f)]

    seqNum = ""
    # gzCheck = glob.glob("%s\w+.fastq.gz" % name)

    # print gzCheck
    if not gzCheck:
    #     # for strain in sampleNames:
        newPath = "%s/%s" % (path, name)
        folderCheck = [f for f in os.listdir(newPath) if re.search("%s\w+.fastq.gz" % name, f)]

    #     # os.chdir("%s/%s" % (path, strain))
    #     # gzFiles = glob.glob("%s*.gz" % name)
    #     # if not gzFiles:
    #     #     pass
        if folderCheck:
            seqNum = re.search("%s_(\w+)_R\d_001.fastq.gz" % name, folderCheck[0])
            moveExtract(name, sorted(folderCheck), path, newPath, seqNum)
    else:
        seqNum = re.search("%s_(\w+)_R\d_001.fastq.gz" % name, gzCheck[0])
        # print name, seqNum.group(1)
    # #     print("Moving and extracting fastq files.")
    # #     for strain in sampleNames:
    # #         os.chdir(path)
    #     # Make the required folders (if necessary)
        make_path("%s/%s" % (path, name))
    # #         # Get the .gz files into a list
    # #         gzFiles = glob.glob("%s_*" % strain)
    # #         if not gzFiles:
    # #             os.chdir("%s/%s" % (path, strain))
    # #             gzFiles = glob.glob("%s/*.gz" % strain)
    # #             if not gzFiles:
    # #                 pass
    # #             else:
        moveExtract(name, sorted(gzCheck), path, path, seqNum)
    #         else:
    #             moveExtract(strain, sorted(gzFiles), path)


def functionsGoNOW(sampleNames, path):
    """Run the functions"""
    unzipping()
    foldererPrepProcesses(sampleNames, path)