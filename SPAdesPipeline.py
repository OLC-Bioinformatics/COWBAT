__author__ = 'akoziol'

# OS commands
import os
# Time module
import time
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

# TODO Move the rMLST and the referenceGenomes folders to a central location
# TODO Figure out how to avoid using absolute paths for called scripts
# TODO Think about getting this pipeline into docker

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/home/blais/PycharmProjects/SPAdesPipeline/2014-09-19")
path = os.getcwd()

# Start time
start = time.time()

# Welcome message
print("Welcome to the CFIA SPAdes Assembly Pipeline.")

# Count the CPUs in the system
cpus = os.popen("awk '/^processor/ { N++} END { print N }' /proc/cpuinfo").read().rstrip()
print("There are %s CPUs in your system" % cpus)


def pipeline():
    """All the functions for running the pipeline"""
    # Import the metadata gathered from GenerateFASTQRunStatistics.xml, RunInfo.xml, and SampleSheet.csv
    print("Extracting metadata from sequencing run.")
    runMetadata, sampleNames, experimentDate = runMetadataOptater.functionsGoNOW()
    # Pre-process archives
    fileExtractionProcessing.functionsGoNOW(sampleNames, path)
    # quakify
    correctedFiles, runTrimMetadata = quakeR.functionsGoNOW(sampleNames, path, runMetadata)
    # SPAdesify
    spadesGoUpper.functionsGoNOW(correctedFiles, path)
    # Typing
    runTrimMLSTMetadata = rMLST_typer.functionsGoNOW(correctedFiles, path, experimentDate, runTrimMetadata)
    # Quasting
    runTrimMLSTAssemblyMetadata = quastR.functionsGoNOW(correctedFiles, path, runTrimMLSTMetadata)
    # Library size estimation
    # print json.dumps(runTrimMLSTMetadata, sort_keys=True, indent=4)
    runTrimMLSTAssemblyInsertMetadata = lse.functionsGoNOW(correctedFiles, path, runTrimMLSTAssemblyMetadata)
    reportR.functionsGoNOW(correctedFiles, runTrimMLSTAssemblyInsertMetadata, path)


# Run the pipeline
pipeline()

print "Elapsed Time: %s seconds" % (time.time() - start)