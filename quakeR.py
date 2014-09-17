__author__ = 'akoziol'

import os
from multiprocessing import Pool
import subprocess
import sys
import re

# Initialise variables
corrected = []


def quakePrepProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakePrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakePool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakePrepArgs.append((name, path))
        # This map function allows for multi-processing
        createQuakePool.map(runQuake, quakePrepArgs)


def runQuake((name, path)):
    """Performs necessary checks and runs quake"""
    # Set up variables to keep commands clean looking
    forward = name + "_R1_001.fastq"
    reverse = name + "_R2_001.fastq"
    statsFile = name + "_R1_001.stats.txt"
    newPath = path + "/" + name
    # Check for the existence of the stats file - hopefully this will be created at the end
    # of the error correction processes
    if not os.path.isfile("%s/%s" % (newPath, statsFile)):
        # Create the list of fastq files required by quake
        fastqList = open("%s/%s_fastqFiles.txt" % (newPath, name), "wb")
        fastqList.write("%s/%s\t%s/%s" % (newPath, forward, newPath, reverse))
        fastqList.close()
        # Quake run command - using a kmer size of 15 because that is what the developers recommended for
        # microbial genomes
        # This is using a hard-coded path, as for some reason, when run within pycharm, quake.py could not
        # be located. Maybe the $PATH needs to be updated?
        # quakeRun = "quake.py -f %s/%s_fastqFiles.txt -k 15 -p 2" % (newPath, name)
        quakeRun = "/home/blais/Bioinformatics/Quake/bin/quake.py -f %s/%s_fastqFiles.txt -k 15 -p 2" % (newPath, name)
        # Run the command

        # subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def completionist(sampleNames, path, runMetadata):
    """This function checks to see which files could not be processed by Quake, and populates a list as appropriate.
    Additionally, it populates a dictionary with the metadata values for the correction"""
    for name in sampleNames:
        statsFileForward = name + "_R1_001.stats.txt"
        statsFileReverse = name + "_R2_001.stats.txt"
        newPath = path + "/" + name
        # Check for the existence of the stats file - hopefully this will be created at the end
        # of the error correction processes
        if not os.path.isfile("%s/%s" % (newPath, statsFileForward)):
            print(name)
            # Populate the dictionary with "N/A" values, as error correction did not occur
            runMetadata[name]["ForwardValidatedReads"] = "N/A"
            runMetadata[name]["ForwardValidatedReads"] = "N/A"
            runMetadata[name]["ForwardCorrectedReads"] = "N/A"
            runMetadata[name]["ForwardTrimmedReads"] = "N/A"
            runMetadata[name]["ForwardTrimmedOnlyReads"] = "N/A"
            runMetadata[name]["ForwardRemovedReads"] = "N/A"
            #
            runMetadata[name]["ReverseValidatedReads"] = "N/A"
            runMetadata[name]["ReverseCorrectedReads"] = "N/A"
            runMetadata[name]["ReverseTrimmedReads"] = "N/A"
            runMetadata[name]["ReverseTrimmedOnlyReads"] = "N/A"
            runMetadata[name]["ReverseRemovedReads"] = "N/A"
        else:
            # Add the strainName to a list of strains that were processed by quake
            corrected.append(name)
            # Open the .stats.txt file with information in the number of Validated, Corrected,
            # Trimmed, Trimmed only, and Removed reads following error correction
            forward = open("%s/%s" % (newPath, statsFileForward), "r")
            for line in forward:
                # Split on the colon
                subline = line.split(": ")
                # Populate the dictionary with the appropriate metadata
                if re.search("Validated", line):
                    runMetadata[name]["ForwardValidatedReads"] = subline[1].strip()
                elif re.search("Corrected", line):
                    runMetadata[name]["ForwardCorrectedReads"] = subline[1].strip()
                elif re.search("Trimmed:", line):
                    runMetadata[name]["ForwardTrimmedReads"] = subline[1].strip()
                elif re.search("only", line):
                    runMetadata[name]["ForwardTrimmedOnlyReads"] = subline[1].strip()
                elif re.search("Removed", line):
                    runMetadata[name]["ForwardRemovedReads"] = subline[1].strip()
            forward.close()
            # Same as with the forward read file
            reverse = open("%s/%s" % (newPath, statsFileReverse), "r")
            for line in reverse:
                subline = line.split(": ")
                if re.search("Validated", line):
                    runMetadata[name]["ReverseValidatedReads"] = subline[1].strip()
                elif re.search("Corrected", line):
                    runMetadata[name]["ReverseCorrectedReads"] = subline[1].strip()
                elif re.search("Trimmed:", line):
                    runMetadata[name]["ReverseTrimmedReads"] = subline[1].strip()
                elif re.search("only", line):
                    runMetadata[name]["ReverseTrimmedOnlyReads"] = subline[1].strip()
                elif re.search("Removed", line):
                    runMetadata[name]["ReverseRemovedReads"] = subline[1].strip()
            reverse.close()
    return runMetadata


def functionsGoNOW(sampleNames, path, runMetadata):
    """Run the functions"""
    print('Performing error correction on fastq files.')
    quakePrepProcesses(sampleNames, path)
    # I don't know why, but it seems that some files don't get processed the first time
    print("\nSecond error correction pass.")
    quakePrepProcesses(sampleNames, path)
    print("\nThese files could not be corrected, and will not be assembled.")
    # Run completionist to determine unprocessable files, and acquire metadata
    runTrimMetadata = completionist(sampleNames, path, runMetadata)
    # Return important variables
    return corrected, runTrimMetadata
