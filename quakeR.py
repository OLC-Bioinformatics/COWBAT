__author__ = 'akoziol'

import os
import subprocess
import sys
import re
import glob

# Initialise variables
corrected = []


def runQuake(sampleName, path):
    """Performs necessary checks and runs quake"""
    for name in sampleName:
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
            quakeRun = "/home/blais/Bioinformatics/Quake/bin/quake.py -f %s/%s_fastqFiles.txt -k 15 -p 24" % (newPath, name)
            # Run the command

            subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')


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
            print(name, " could not be corrected, and will not be assembled.")
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
    # Removed the multiprocessing aspect of this function - it seemed to be unreliable.
    # Sometimes, fastq files with more data would not be corrected.
    runQuake(sampleNames, path)
    # Run completionist to determine unprocessable files, and acquire metadata
    runTrimMetadata, correctedList = completionist(sampleNames, path, runMetadata)
    # Clean up tmp files
    tmpFileRemover(path)
    # Return important variables
    return correctedList, runTrimMetadata
