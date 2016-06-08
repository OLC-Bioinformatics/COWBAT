__author__ = 'akoziol'

import os
import subprocess
import sys
import re
import glob
from multiprocessing import Pool
import json
import metadataFiller
import jsonReportR
# Initialise variables
corrected = []

def quakePrepProcesses(sampleName, path, fLength, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakePrepArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakePool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakePrepArgs.append((name, path, fLength, metadata, commands))
        # This map function allows for multi-processing
        output = createQuakePool.map(runQuakeMP, quakePrepArgs)
    return output


def runQuakeMP((name, path, fLength, metadata, commands)):
    """Multiprocessed version of runQuake"""
    # Initialise the variables for the forward and reverse reads
    forward = name + "_R1_001.fastq"
    reverse = name + "_R2_001.fastq"
    newPath = path + "/" + name
    # This glob is used to allow for the processing of files that do not match the established naming conventions
    forwardGlob = glob.glob("%s/*1.fastq" % newPath)
    reverseGlob = glob.glob("%s/*2.fastq" % newPath)
    # The counts file and its associated stats are required
    countsFile = name + "_counts.txt"
    countSize = 0
    if os.path.isfile("%s/%s" % (newPath, countsFile)):
        countSize = os.stat("%s/%s" % (newPath, countsFile)).st_size
    # Create the list of fastq files required by quake
    fastqList = open("%s/%s_fastqFiles.txt" % (newPath, name), "wb")
    # fastqList.write("%s\t%s" % (forward, reverse))
    if fLength > 50:
        fastqList.write("%s/%s\t%s/%s" % (newPath, forward, newPath, reverse))
    else:
        fastqList.write("%s/%s" % (newPath, reverse))
    fastqList.close()
    corFile = glob.glob("%s/*cor.fastq" % newPath)
    # if not os.path.isfile("%s/%s" % (newPath, countsFile))
    if not commands[name]["QuakeQmersCorrectCommand"] and countSize == 0 and not corFile:
        # if not os.path.isfile("%s/%s" % (newPath, countsFile)) or countSize == 0 and len(reverseGlob) > 0:
            # Looks for the reverse file, if it doesn't exist, try again, by a more general regex glob
        if not os.path.isfile("%s/%s" % (newPath, reverse)) and forwardGlob and reverseGlob:
            forwardPath = forwardGlob[0]
            forward = os.path.split(forwardPath)[1]
            reversePath = reverseGlob[0]
            reverse = os.path.split(reversePath)[1]
        # Quake run command - using a kmer size of 15 because that is what the developers recommended for
        # microbial genomes. Also using 24 processors.
        #   2>/dev/null -q 33 -k 15
        if fLength > 50:
            quakeRun = "cat %s/%s %s/%s | count-qmers -q 33 -k 15 > %s/%s" \
                       % (newPath, forward, newPath, reverse, newPath, countsFile)
            metadata[name]["7.PipelineCommands"]["QuakeQmersCorrectCommand"] = quakeRun
        else:
            quakeRun = "cat %s/%s | count-qmers -q 33 -k 15 > %s/%s" \
                       % (newPath, reverse, newPath, countsFile)
            metadata[name]["7.PipelineCommands"]["QuakeQmersCorrectCommand"] = quakeRun
        # Run the command
        os.system(quakeRun)
        sys.stdout.write('.')
        return metadata
        # else:
        #     sys.stdout.write('.')
        #     metadata[name]["7.PipelineCommandsCommands"]["QuakeQmersCommand"] = commands[name]["QuakeQmersCorrectCommand"]
        #     return metadata
    else:
        sys.stdout.write('.')
        metadata[name]["7.PipelineCommands"]["QuakeQmersCorrectCommand"] = commands[name]["QuakeQmersCorrectCommand"]
        return metadata


def quakeCutOffProcesses(sampleName, path, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakeCutArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakeCutPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakeCutArgs.append((name, path, metadata, commands))
        # This map function allows for multi-processing
        output = createQuakeCutPool.map(cutQuakeMP, quakeCutArgs)
    return output


def cutQuakeMP((name, path, metadata, commands)):
    """Multiprocessed version of runQuake"""
    countSize = 0
    countsFile = name + "_counts.txt"
    newPath = path + "/" + name
    if os.path.isfile("%s/%s" % (newPath, countsFile)):
        countSize = os.stat("%s/%s" % (newPath, countsFile)).st_size
    corFile = glob.glob("%s/*cor.fastq" % newPath)
    # Check for the existence of the assembled contigs - if this file exists, then skip
    if not commands[name]["QuakeCutCommand"] and countSize != 0 and not corFile:
        # If cutoff.txt doesn't exist, or counts.txt was not populated properly (or yet), then
        # if not os.path.isfile("%s/cutoff.txt" % newPath) or countSize == 0:
        # I made edits to the cov_model.py script, and the R script called by cov_model.py to allow
        # for multi-processing. Essentially, I modified the scripts to include a path variable, so that
        # all files searched for and created are in 'newPath' instead of the path
        quakeCut = "cov_model.py --path %s %s/%s 2>/dev/null" % (newPath, newPath, countsFile)
        # print quakeCut
        # Run the command
        os.system(quakeCut)
        metadata[name]["7.PipelineCommands"]["QuakeCutCommand"] = quakeCut
        sys.stdout.write('.')
        return metadata
        # else:
        #     sys.stdout.write('.')
        #     metadata[name]["7.PipelineCommands"]["QuakeCutCommand"] = commands[name]["QuakeCutCommand"]
        #     return metadata
    else:
        metadata[name]["7.PipelineCommands"]["QuakeCutCommand"] = commands[name]["QuakeCutCommand"]
        sys.stdout.write('.')
        return metadata


def quakeCorrectProcesses(sampleName, path, fLength, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    quakeCorrectArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quakeR':
        createQuakeCorrectPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            quakeCorrectArgs.append((name, path, fLength, metadata, commands))
        # This map function allows for multi-processing
        output = createQuakeCorrectPool.map(correctQuakeMP, quakeCorrectArgs)
    # print json.dumps(output, sort_keys=True, indent=4, separators=(',', ': '))
    return output


def correctQuakeMP((name, path, fLength, metadata, commands)):
    """Multiprocessed version of runQuake"""
    newPath = path + "/" + name
    cutoffSize = 0
    cutoffFile = "%s/cutoff.txt" % newPath
    if os.path.isfile(cutoffFile):
        cutoffSize = os.stat(cutoffFile).st_size

    corFile = glob.glob("%s/*cor.fastq" % newPath)
    # os.path.isfile("%s/cutoff.txt" % newPath)
    # and cutoffSize != 0
    if not commands[name]["QuakeCorrectCommand"] and not corFile:
        # As part of the assembly of GeneSippr data, only one of the two paired end reads are necessary
        if fLength > 50:
            necessaryNoCorFiles = 2
        else:
            necessaryNoCorFiles = 1
        # if os.path.isfile(cutoffFile) and cutoffSize == 0 and len(corFile) < necessaryNoCorFiles:
        if os.path.isfile('%s/cutoff.txt' % newPath):
            cutoff = open('%s/cutoff.txt' % newPath).readline().rstrip()
            if cutoff == "Inf":
                cutoff = 1
        else:
            # Otherwise set the cutoff to 1 if the cutoff command failed
            cutoff = 1
        # microbial genomes. Also using 24 processors.
        #  2>/dev/null  -q 33 -k 15
        quakeCorrect = "correct -f %s/%s_fastqFiles.txt -k 15 -m %s/%s_counts.txt -c %s -p 24" % (newPath, name, newPath, name, cutoff)
        print quakeCorrect
        # quakeRun = "quake.py -f %s/%s_fastqFiles.txt -k 15 -p 24" % (newPath, name)
        # Run the command
        # subprocess.call(quakeRun, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # print quakeCorrect
        os.system(quakeCorrect)
        metadata[name]["7.PipelineCommands"]["QuakeCorrectCommand"] = quakeCorrect
        sys.stdout.write('.')
        return metadata
        # else:
        #     metadata[name]["7.PipelineCommands"]["QuakeCorrectCommand"] = commands[name]["QuakeCorrectCommand"]
        #     sys.stdout.write('.')
        #     return metadata
    else:
        metadata[name]["7.PipelineCommands"]["QuakeCorrectCommand"] = commands[name]["QuakeCorrectCommand"]
        sys.stdout.write('.')
        return metadata


def completionist(sampleNames, path, runMetadata, fLength):
    """This function checks to see which files could not be processed by Quake, and populates a list as appropriate.
    Additionally, it populates a dictionary with the metadata values for the correction"""
    for name in sampleNames:
        # Populate the number of reads
        newPath = path + "/" + name
        runMetadata[name]["4.Correction"]["kmerSize"] = 15
        statsFileForward = name + "_R1_001.stats.txt"
        statsFileReverse = name + "_R2_001.stats.txt"
        if not os.path.isfile("%s/%s" % (newPath, statsFileReverse)):
            statsGlob = glob.glob("%s/*.stats.txt" % newPath)
            if statsGlob:
                if fLength > 50:
                    statsFileForwardPath = sorted(statsGlob)[0]
                    statsFileForward = os.path.split(statsFileForwardPath)[1]
                statsFileReversePath = sorted(statsGlob)[1]
                statsFileReverse = os.path.split(statsFileReversePath)[1]
        # Check for the existence of the stats file - hopefully this will be created at the end
        # of the error correction processes
        if not os.path.isfile("%s/%s" % (newPath, statsFileReverse)):
            print "\n%s could not be corrected, and will not be assembled." % name
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
            if fLength > 50:
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
            else:
                runMetadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
                runMetadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
                runMetadata[name]["4.Correction"]["ForwardCorrectedReads"] = "N/A"
                runMetadata[name]["4.Correction"]["ForwardTrimmedReads"] = "N/A"
                runMetadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"] = "N/A"
                runMetadata[name]["4.Correction"]["ForwardRemovedReads"] = "N/A"
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


def tmpFileRemover(path, correctedList):
    """Removes temporary files cluttering up the path"""
    if os.path.isfile("%s/r.log" % path):
        os.remove("%s/r.log" % path)
    tmpFiles = glob.glob("%s/*.txt*" % path)
    for file in tmpFiles:
        if not re.search("indexingQC.txt", file):
            os.remove(file)
    for name in correctedList:
        newPath = path + "/" + name
        fastq = glob.glob("%s/*_001.fastq" % newPath)
        for files in fastq:
            if os.path.isfile(files):
                os.remove(files)
        counts = "%s/%s_counts.txt" % (newPath, name)
        if os.path.isfile(counts):
            os.remove(counts)
        qcts = "%s/%s_fastqFiles.txt.qcts" % (newPath, name)
        if os.path.isfile(qcts):
            os.remove(qcts)


def functionsGoNOW(sampleNames, path, runMetadata, fLength, commands):
    """Run the functions"""
    print('\nPerforming error correction on fastq files.')
    # Removed the multiprocessing aspect of this function - it seemed to be unreliable.
    # Sometimes, fastq files with more data would not be corrected.
    os.chdir(path)
    print "Preparing fastq files for processing"
    prepList = quakePrepProcesses(sampleNames, path, fLength, runMetadata, commands)
    prepMetadata = metadataFiller.filler(runMetadata, prepList)
    print "Determining cut-off values for error correction"
    cutoffList = quakeCutOffProcesses(sampleNames, path, prepMetadata, commands)
    cutoffMetadata = metadataFiller.filler(prepMetadata, cutoffList)
    print "Correcting errors"
    correctList = quakeCorrectProcesses(sampleNames, path, fLength, cutoffMetadata, commands)
    correctMetadata = metadataFiller.filler(runMetadata, correctList)
    # runQuake(sampleNames, path)
    os.chdir(path)
    # Run completionist to determine unprocessable files, and acquire metadata
    runTrimMetadata, correctedList = completionist(sampleNames, path, correctMetadata, fLength)
    # Clean up tmp files
    tmpFileRemover(path, correctedList)
    # Return important variables
    jsonReportR.jsonR(correctedList, path, runTrimMetadata, "Collection")
    return correctedList, runTrimMetadata
