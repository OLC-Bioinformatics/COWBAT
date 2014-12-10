__author__ = 'akoziol'

import json
import errno
import os
import shutil
import re
from types import *
import glob

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def functionsGoNOW(sampleNames, metadata, path):
    """As outputting a JSON file is so straightforward, helper functions were unnecessary"""
    print "\nCreating reports."
    reportPath = "%s/reports" % path
    make_path(reportPath)
    combinedReport = open("%s/CombinedMetadataReport.csv" % reportPath, "wb")
    headings = ["fileName", "SampleName", "N50", "NumContigs", "TotalLength", "MeanInsertSize", "averageDepthofCov",
                "referenceGenome", "NumIdenticalAlleles", "rMLSTSequenceType", "rMLSTIdenticalAlleles", "Date",
                "percentGC", "Investigator", "TotalClustersinRun", "LengthofFirstRead", "LengthofSecondRead", "Project"]
    reportHeadings = ["1.General", "2.Assembly", "3.Run", "4.Correction", "5.rMLST", "6.rMLSTmatchestoRef"]
    combinedReport.write("\t".join(headings))
    combinedReport.write("\n")
    for name in sampleNames:
        newPath = path + "/" + name
        reportName = "%s_metadataReport.json" % name
        JSONreport = open("%s/%s" % (newPath, reportName), "wb")
        # Print the JSON data to file
        output = json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
        JSONreport.write(output)
        JSONreport.close()
        # Move all the reports to a common directory
        shutil.copy("%s/%s" % (newPath, reportName), reportPath)
        #
        for heading in headings:
            value = "NA"
            for rHeading in reportHeadings:
                if type(metadata[name][rHeading][heading]) is StringType or type(metadata[name][rHeading][heading]) \
                        is FloatType or type(metadata[name][rHeading][heading]) is IntType:
                    # print name, heading, rHeading, metadata[name][rHeading][heading], type(metadata[name][rHeading][heading])
                    value = str(metadata[name][rHeading][heading])
                    # combinedReport.write(str(metadata[name][rHeading][heading]))
                    # combinedReport.write("\t")
                # else:

            combinedReport.write(value)
            combinedReport.write("\t")
        combinedReport.write("\n")

        # file clean-up
        dbm = "%s/%s_fastqFiles.dbm" % (newPath, name)
        if os.path.isfile(dbm):
            os.remove(dbm)
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