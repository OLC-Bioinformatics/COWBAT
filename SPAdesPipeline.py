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
# Custom Script to perform quake error corrections on the reads
import quakeR
# Run SPAdes
import spadesGoUpper
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

# TODO Move the rMLST and the referenceGenomes folders to a central location
# TODO Figure out how to avoid using absolute paths for called scripts
# TODO Locate MiSeq mount and copy files to a new folder in the nas using the date of the run as the file name
# TODO Think about getting this pipeline into docker
# TODO Look at filtering at 10X coverage in SpadesGoUpper
# TODO Add the ability to dynamically determine the versions of all third-party software

refFilePath = "/spades_pipeline/SPAdesPipelineFiles"

# Start time
start = time.time()

# Welcome message
print("Welcome to the CFIA SPAdes Assembly Pipeline.")

# Count the CPUs in the system
cpus = os.popen("awk '/^processor/ { N++} END { print N }' /proc/cpuinfo").read().rstrip()
print("There are %s CPUs in your system" % cpus)

def pipeline():
    """All the functions for running the pipeline"""
    # Off-hours module
    # offHours.run()
    path = os.getcwd()
    # path = "/media/nas/akoziol/WGS_Spades/2014-05-30"
    # path = "/media/nas/akoziol/WGS_Spades/2015-02-10"
    # os.chdir(path)
    # Import the metadata gathered from GenerateFASTQRunStatistics.xml, RunInfo.xml, and SampleSheet.csv
    print("Extracting metadata from sequencing run.")
    runMetadata, sampleNames, experimentDate, fLength = runMetadataOptater.functionsGoNOW(path)
    # print json.dumps(runMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Pre-process archives
    fileExtractionProcessing.functionsGoNOW(sampleNames, path)
    # quakify
    correctedFiles, runTrimMetadata = quakeR.functionsGoNOW(sampleNames, path, runMetadata, fLength)
    # print json.dumps(runTrimMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # SPAdesify
    runTrimAssemblyMetadata, assembledFiles = spadesGoUpper.functionsGoNOW(correctedFiles, path, runTrimMetadata, fLength)
    # print json.dumps(runTrimAssemblyMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Typing
    runTrimAssemblyMLSTMetadata = rMLST_typer.functionsGoNOW(assembledFiles, path, experimentDate, runTrimAssemblyMetadata, refFilePath)
    # print json.dumps(runTrimAssemblyMLSTMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Quasting
    # runTrimMLSTMetadata
    runTrimAssemblyMLSTQuastMetadata = quastR.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTMetadata)
    # print json.dumps(runTrimAssemblyMLSTQuastMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Library size estimation
    runTrimAssemblyMLSTQuastInsertMetadata = lse.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTQuastMetadata)
    # Univec screening

    # Mobile element screening
    # print json.dumps(runTrimMLSTAssemblyInsertMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # GeneSeeking
    runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata = geneSeekr.functionsGoNOW(assembledFiles, path, runTrimAssemblyMLSTQuastInsertMetadata, refFilePath)
    # print json.dumps(runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    # Generate the final metadata reports
    reportR.functionsGoNOW(assembledFiles, runTrimAssemblyMLSTQuastInsertgeneSeekrMetadata, path)


# Run the pipeline
pipeline()

print "\nElapsed Time: %s seconds" % (time.time() - start)