__author__ = 'akoziol'

import os
import subprocess
import sys
import re
import glob
from multiprocessing import Pool
import json

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
        createQuakePool.map(runQuakeMP, quakePrepArgs)


def runQuakeMP((name, path)):
    """Multiprocessed version of runQuake"""
    forward = name + "_R1_001.fastq"
    reverse = name + "_R2_001.fastq"
    countsFile = name + "_counts.txt"
    newPath = path + "/" + name
    # Check for the existence of the stats file - hopefully this will be created at the end
    # of the error correction processes
    # os.remove("%s/%s_fastqFiles.txt" % (newPath, name))
    # print "%s/%s_fastqFiles.txt" % (newPath, name)
    # print newPath, statsFile
    if not os.path.isfile("%s/%s" % (newPath, countsFile)) or not os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, name)):
        # Create the list of fastq files required by quake
        fastqList = open("%s/%s_fastqFiles.txt" % (newPath, name), "wb")
        # fastqList.write("%s\t%s" % (forward, reverse))
        fastqList.write("%s/%s\t%s/%s" % (newPath, forward, newPath, reverse))
        fastqList.close()
        # Quake run command - using a kmer size of 15 because that is what the developers recommended for
        # microbial genomes. Also using 24 processors.
        quakeRun = "cat %s/%s %s/%s | count-qmers -q 33 -k 15 > %s/%s" % (newPath, forward, newPath, reverse, newPath, countsFile)
        # quakeRun = "quake.py -f %s/%s_fastqFiles.txt -k 15 -p 24" % (newPath, name)
        # Run the command
        # subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        os.system(quakeRun)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def quakeCutOffProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakeCutArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakeCutPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakeCutArgs.append((name, path))
        # This map function allows for multi-processing
        createQuakeCutPool.map(cutQuakeMP, quakeCutArgs)


def cutQuakeMP((name, path)):
    """Multiprocessed version of runQuake"""
    countsFile = name + "_counts.txt"
    newPath = path + "/" + name
    # Check for the existence of the stats file - hopefully this will be created at the end
    # of the error correction processes
    # os.remove("%s/%s_fastqFiles.txt" % (newPath, name))
    # print "%s/%s_fastqFiles.txt" % (newPath, name)
    # print newPath, statsFile
    if not os.path.isfile("%s/cutoff.txt" % newPath):
        # I made edits to the cov_model.py script, and the R script called by cov_model.py to allow
        # for multi-processing. Essentially, I modified the scripts to include a path variable, so that
        # all files searched for and created are in 'newPath' instead of the path
        quakeCut = "cov_model.py --path %s %s/%s" % (newPath, newPath, countsFile)
        # Run the command
        # subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        os.system(quakeCut)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def quakeCorrectProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakeCorrectArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakeCorrectPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakeCorrectArgs.append((name, path))
        # This map function allows for multi-processing
        createQuakeCorrectPool.map(correctQuakeMP, quakeCorrectArgs)


def correctQuakeMP((name, path)):
    """Multiprocessed version of runQuake"""
    newPath = path + "/" + name
    # Check for the existence of the stats file - hopefully this will be created at the end
    # of the error correction processes
    # os.remove("%s/%s_fastqFiles.txt" % (newPath, name))
    # print "%s/%s_fastqFiles.txt" % (newPath, name)
    # print newPath, statsFile
    if os.path.isfile("%s/cutoff.txt" % newPath) and not os.path.isfile("%s/%s_R1_001.stats.txt" % (newPath, name)):
        cutoff = open('%s/cutoff.txt' % newPath).readline().rstrip()
        # microbial genomes. Also using 24 processors.
        quakeCorrect = "correct -f %s/%s_fastqFiles.txt -k 15 -q 33 -m %s/%s_counts.txt -c %s -p 24" % (newPath, name, newPath, name, cutoff)
        # quakeRun = "quake.py -f %s/%s_fastqFiles.txt -k 15 -p 24" % (newPath, name)
        # Run the command
        # subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        os.system(quakeCorrect)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')

# def runQuake(sampleName, path):
#     """Performs necessary checks and runs quake"""
#     os.chdir(path)
#     for name in sampleName:
#         # Set up variables to keep commands clean looking
#         forward = name + "_R1_001.fastq"
#         reverse = name + "_R2_001.fastq"
#         statsFile = name + "_R1_001.stats.txt"
#         newPath = path + "/" + name
#         os.chdir(newPath)
#         # os.remove("%s_fastqFiles.txt" % name)
#         # Check for the existence of the stats file - hopefully this will be created at the end
#         # of the error correction processes
#         if not os.path.isfile(statsFile):
#             # Create the list of fastq files required by quake
#             fastqList = open("%s_fastqFiles.txt" % name, "wb")
#             fastqList.write("%s\t%s" % (forward, reverse))
#             fastqList.close()
#             # Quake run command - using a kmer size of 15 because that is what the developers recommended for
#             # microbial genomes. Also using 24 processors.
#             quakeRun = "quake.py -f %s_fastqFiles.txt -k 15 -p 24" % name
#             # Run the command
#             subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
#             sys.stdout.write('.')
#         else:
#             sys.stdout.write('.')


def completionist(sampleNames, path, runMetadata):
    """This function checks to see which files could not be processed by Quake, and populates a list as appropriate.
    Additionally, it populates a dictionary with the metadata values for the correction"""
    for name in sampleNames:
        # Populate the number of reads
        statsFileForward = name + "_R1_001.stats.txt"
        statsFileReverse = name + "_R2_001.stats.txt"
        newPath = path + "/" + name
        # Check for the existence of the stats file - hopefully this will be created at the end
        # of the error correction processes
        if not os.path.isfile("%s/%s" % (newPath, statsFileForward)):
            print "\n%s could not be corrected, and will not be assembled." % name
            # corrected.append(name)
            # Populate the dictionary with "N/A" values, as error correction did not occur
            runMetadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ForwardCorrectedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ForwardTrimmedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ForwardRemovedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ReverseValidatedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ReverseCorrectedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ReverseTrimmedReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ReverseTrimmedOnlyReads"] = "N/A"
            runMetadata[name]["4.Correction"]["ReverseRemovedReads"] = "N/A"
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
                    runMetadata[name]["4.Correction"]["ForwardValidatedReads"] = subline[1].strip()
                elif re.search("Corrected", line):
                    runMetadata[name]["4.Correction"]["ForwardCorrectedReads"] = subline[1].strip()
                elif re.search("Trimmed:", line):
                    runMetadata[name]["4.Correction"]["ForwardTrimmedReads"] = subline[1].strip()
                elif re.search("only", line):
                    runMetadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"] = subline[1].strip()
                elif re.search("Removed", line):
                    runMetadata[name]["4.Correction"]["ForwardRemovedReads"] = subline[1].strip()
            forward.close()
            # Same as with the forward read file
            reverse = open("%s/%s" % (newPath, statsFileReverse), "r")
            for line in reverse:
                subline = line.split(": ")
                if re.search("Validated", line):
                    runMetadata[name]["4.Correction"]["ReverseValidatedReads"] = subline[1].strip()
                elif re.search("Corrected", line):
                    runMetadata[name]["4.Correction"]["ReverseCorrectedReads"] = subline[1].strip()
                elif re.search("Trimmed:", line):
                    runMetadata[name]["4.Correction"]["ReverseTrimmedReads"] = subline[1].strip()
                elif re.search("only", line):
                    runMetadata[name]["4.Correction"]["ReverseTrimmedOnlyReads"] = subline[1].strip()
                elif re.search("Removed", line):
                    runMetadata[name]["4.Correction"]["ReverseRemovedReads"] = subline[1].strip()
            reverse.close()
    return runMetadata, corrected


def tmpFileRemover(path):
    """Removes temporary files cluttering up the path"""
    if os.path.isfile("%s/r.log" % path):
        os.remove("%s/r.log" % path)
    tmpFiles = glob.glob("%s/*.txt*" % path)
    for file in tmpFiles:
        os.remove(file)


def functionsGoNOW(sampleNames, path, runMetadata):
    """Run the functions"""
    print('\nPerforming error correction on fastq files.')
    # print json.dumps(runMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Removed the multiprocessing aspect of this function - it seemed to be unreliable.
    # Sometimes, fastq files with more data would not be corrected.
    os.chdir(path)
    quakePrepProcesses(sampleNames, path)
    quakeCutOffProcesses(sampleNames, path)
    quakeCorrectProcesses(sampleNames, path)
    # runQuake(sampleNames, path)
    os.chdir(path)
    # Run completionist to determine unprocessable files, and acquire metadata
    runTrimMetadata, correctedList = completionist(sampleNames, path, runMetadata)
    # Clean up tmp files
    tmpFileRemover(path)
    # Return important variables
    return correctedList, runTrimMetadata
