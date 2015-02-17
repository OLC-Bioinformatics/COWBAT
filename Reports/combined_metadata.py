#! /usr/bin/env python
__author__ = 'akoziol'

import json
import glob
import os
import sys
import re

metadataResults = ["fileName", "MedianInsertSize", "InsertSizeStDev", "rMLSTSequenceType", "NumIdenticalAlleles",
                   "referenceGenome", "rMLSTmatchestoRef", "N50", "LargestContig", "NumContigs", "TotalLength"]

path = os.getcwd()

files = glob.glob("*.json")

# jsonFile = json.loads(open(files[0]).read())
jsonReadFile = open(files[0], "r").readlines()

# print json.dumps(jsonFile["1.General"]["InsertSizeStDev"]).replace('"', '')
# print jsonFile["1.General"]["rMLSTSequenceType"]["sequenceType"]
# for key in jsonFile["1.General"]["rMLSTSequenceType"]["sequenceType"]:
    # print key
for metadata in metadataResults:
    for line in jsonReadFile:
        if re.search(metadata, line):
            print metadata, line.split(":")[1].replace(",", "").replace('"', '').strip()