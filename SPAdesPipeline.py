#! /usr/bin/env python
__author__ = 'akoziol'

# OS commands
import os
# Time module
import time
# Custom script to pull the data from the MiSeq
import offHours
# Custom script for pulling metadata from sequencing run reports/files
import runMetadataOptater
# Custom script for moving and/or extracting archived files
import fileExtractionProcessing
import singleReadsExtraction
# Custom Script to perform quake error corrections on the reads
import quakeR
# Run SPAdes
import spadesGoUpper
import spadesUpGoer
# Perform typing using rMLST to determine best reference genome
import rMLST_typer
# quastR
import quastR
# Library size estimation
import lse
# Create a JSON report
import reportR
# Perform GeneSeeking analysis
import geneSeekr
import json
import commandFile
# Argument parser for user-inputted values, and a nifty help menu
from argparse import ArgumentParser

#Parser for arguments
parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
parser.add_argument('-v', '--version', action='version', version='%(prog)s commit 5c63c82cd21e9e1883cd0d69e91a93ccff329b9c')
parser.add_argument('-p', '--path', required=True, help='Specify path')
parser.add_argument('-n', '--numReads', required=False, default=2,
                    help='Specify the number of reads. Paired-reads: 2, unpaired-reads: 1. Default is paired-end (2)')
parser.add_argument('-t', '--threads', required=False, help='Number of threads to use. Defaults to the number of '
                    'cores in your system')
parser.add_argument('-o', '--offHours', required=False, help='Optionally run the off-hours module, which will '
                    'search the MiSeq for folders with data that has not been assembled')
parser.add_argument('-r', '--referenceFilePath', required=False, default="/spades_pipeline/SPAdesPipelineFiles",
                    help="Provide the location of the folder containing the pipeline accessory files (reference genomes, MLST data, etc.")

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
path = args['path']
cpus = args['threads']
numReads = args['numReads']
offhours = args['offHours']
refFilePath = args['referenceFilePath']
os.chdir(path)

# TODO Look at filtering at 10X coverage in SpadesGoUpper
# TODO Save discarded contigs in file

# I haven't implemented this, as sometimes, it's better to have coverage below 10X than deleting all the contigs

# TODO add pxz of fastq files at end

# refFilePath = "/spades_pipeline/SPAdesPipelineFiles"

# Start time
start = time.time()

# Welcome message
print("Welcome to the CFIA SPAdes Assembly Pipeline.")

# Count the CPUs in the system
if not cpus:
    cpus = os.popen("awk '/^processor/ { N++} END { print N }' /proc/cpuinfo").read().rstrip()
print("There are %s CPUs in your system" % cpus)


def pipeline():
    """All the functions for running the pipeline"""
    global path
    # Off-hours module - run if -o flag present
    if offhours:
        offHours.run()
        path = os.getcwd()
    # Import the metadata gathered from GenerateFASTQRunStatistics.xml, RunInfo.xml, and SampleSheet.csv
    print("Extracting metadata from sequencing run.")
    runMetadata, sampleNames, experimentDate, fLength = runMetadataOptater.functionsGoNOW(path)
    # print json.dumps(runMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    commands, versions, reads, commandMetadata = commandFile.commands(sampleNames, path, numReads, runMetadata)
    # Pre-process archives if int(numReads) == 1:
    if int(numReads) == 1 or reads == "single":
        singleReadsExtraction.functionsGoNOW(sampleNames, path)
        # Run the special single read Spades assembly function
        runTrimAssemblyMetadata, assembledFiles = spadesUpGoer.functionsGoNow(sampleNames, path, commandMetadata, fLength, commands)
    else:
        fileExtractionProcessing.functionsGoNOW(sampleNames, path)
        # quakify
        correctedFiles, runTrimMetadata = quakeR.functionsGoNOW(sampleNames, path, commandMetadata, fLength, commands)
        # print json.dumps(runTrimMetadata, sort_keys=True, indent=4, separators=(',', ': '))
        # SPAdesify
        runTrimAssemblyMetadata, assembledFiles = spadesGoUpper.functionsGoNOW(correctedFiles, path, runTrimMetadata, fLength, commands)
    #     # print json.dumps(runTrimAssemblyMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Typing
    runTrimAssemblyMLSTMetadata = rMLST_typer.functionsGoNOW(assembledFiles, path, experimentDate, runTrimAssemblyMetadata, refFilePath)
    # Quasting
    # # runTrimMLSTMetadata
    runTrimAssemblyMLSTQuastMetadata = quastR.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTMetadata, commands)
    # print json.dumps(runTrimAssemblyMLSTQuastMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Library size estimation
    runTrimAssemblyMLSTQuastInsertMetadata = lse.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTQuastMetadata, commands)
    # print json.dumps(runTrimAssemblyMLSTQuastInsertMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Mobile element screening
    # Not implemented yet
    # GeneSeeking
    runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata = geneSeekr.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTQuastInsertMetadata, refFilePath, commands)
    # print json.dumps(runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata["996-Old"]["7.PipelineCommands"], sort_keys=True, indent=4, separators=(',', ': '))
    # # Generate the final metadata reports
    reportR.functionsGoNOW(assembledFiles, runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata, path)


# Run the pipeline
pipeline()

print "\nElapsed Time: %s seconds" % (time.time() - start)