__author__ = 'akoziol'

import re
import os
from collections import defaultdict
import pprint
import json

# Initialise variables
global totalClustersPF
reads = []
samples = []

# global investigator

# Import ElementTree - try first to import the faster C version, if that doesn't
# work, try to import the regular version
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
reportData = defaultdict(make_dict)
sampleData = defaultdict(make_dict)

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/media/nas1/akoziol/Pipeline_development/SPAdesPipelineSandbox")
path = os.getcwd()


def parseSampleSheet():
    """Parses the sample sheet (SampleSheet.csv) to determine certain values
    important for the creation of the assembly report"""
    global investigator
    global experiment
    global date
    global forwardLength
    global reverseLength
    global adapter
    sampleSheet = open("SampleSheet.csv", "r")
    for line in sampleSheet:
        line.strip()
        data = line.split(",")
        if re.search("Investigator", line):
            investigator = data[1].rstrip()
        elif re.search("Experiment", line):
            experiment = data[1].rstrip().replace("  ", " ")
        elif re.search("Date", line):
            date = data[1].rstrip()
        # Here's a Perl-like solution for reading lines after a regex match
        # Perform the regex
        elif re.search("Reads", line):
            # Now read all sublines going forward
            for subline in sampleSheet:
                # Stop reading once "Settings" is matched
                if re.search("Settings", subline):
                    break
                reads.append(subline)
            forwardLength = reads[0].rstrip()
            reverseLength = reads[1].rstrip()
        if re.search("Adapter", line):
            adapter = data[1].rstrip()
        elif re.search("Sample_ID", line):
            for subline in sampleSheet:
                subdata = subline.split(",")
                # Capture      Sample_ID	          Sample_Name	         I7_Index_ID	        index	         I5_Index_ID	          index2	       Sample_Project
                sampleData[subdata[0].rstrip()][subdata[1].rstrip()][subdata[4].rstrip()][subdata[5].rstrip()][subdata[6].rstrip()][subdata[7].rstrip()] = subdata[8].rstrip()
        # elif re.search("Settings", line):
            # print line






def parseRunStats():
    """Parses the XML run statistics (GenerateFASTQRunStatistics.xml)"""
    pass

parseSampleSheet()
print investigator, experiment, date, forwardLength, reverseLength, adapter

print json.dumps(sampleData, sort_keys=True, indent=4)