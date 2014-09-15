__author__ = 'akoziol'

# Regex
import re
# OS commands
import os
# Perl-style dictionaries
from collections import defaultdict
# Prints variables in an easy-to-ready JSON format
import json
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
# Time module
import time
# Custom script for pulling metadata from sequencing run reports/files
import runMetadataOptater
# Custom script for moving and/or extracting archived files
import fileExtractionProcessing

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/media/nas1/akoziol/Pipeline_development/SPAdesPipelineSandbox")
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
    runMetadata, sampleNames = runMetadataOptater.functionsGoNOW()
    # Pre-process archives
    fileExtractionProcessing.functionsGoNOW(sampleNames, path)

# Run the pipeline
pipeline()

print "\nElapsed Time: %s seconds" % (time.time() - start)
