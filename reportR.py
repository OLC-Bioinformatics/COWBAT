__author__ = 'akoziol'

import json
import errno
import os
import shutil
import re
from types import *
import glob
import jsonReportR
import versionIdentifyR

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def reportWriter(sampleNames, metadata, path):
    """As outputting a JSON file is so straightforward, helper functions were unnecessary"""
    print "\nCreating reports."
    reportPath = "%s/reports" % path
    # Grab the name of the analysis from the path variable
    folderName = path.split("/")[-1]
    # Get the path of the parental folder. This is where subfolders containing copies of
    # all the assemblies and reports will be stored
    repositoryPath = os.path.dirname(path)
    # Create the appropriate folders (if necessary) for storing the appropriate files
    make_path("%s/AssemblyData/Assemblies" % repositoryPath)
    make_path("%s/AssemblyData/JsonReports" % repositoryPath)
    make_path("%s/AssemblyData/SummaryReports" % repositoryPath)
    make_path(reportPath)
    combinedReport = open("%s/%s_CombinedMetadataReport.tsv" % (reportPath, folderName), "wb")
    # The headings will be used to search the json file and return only these particular data for
    # inclusion in the metadata spreadsheet
    headings = ["SampleName", "fileName", "N50", "NumContigs", "TotalLength", "MeanInsertSize", "averageDepthofCov",
                "referenceGenome", "NumIdenticalAlleles", "rMLSTSequenceType", "rMLSTIdenticalAlleles", "SequencingDate",
                "NumPredictedGenes", "NumPredictedGenes>500bp", "NumPredictedGenes>1000bp", "NumPredictedGenes>3000bp",
                "geneSeekrProfile", "verotoxinProfile", "MLST_sequenceType", "percentGC", "Investigator",
                "NumberOfClustersPF", "PercentOfClusters", "TotalClustersinRun", "LengthofFirstRead",
                "LengthofSecondRead", "Project", "PipelineVersion"]
    # General headings within the metadata json file
    reportHeadings = ["1.General", "2.Assembly", "3.Run", "4.Correction",
                      "5.rMLST", "6.rMLSTmatchestoRef", "7.PipelineCommands", "8.PipelineVersions"]
    combinedReport.write("\t".join(headings))
    combinedReport.write("\n")
    # Use jsonReportR to fill out the final json metadata report
    jsonReportR.jsonR(sampleNames, path, metadata, "Report")
    # Fill the combined metadata report spreadsheet with the appropriate data
    for name in sampleNames:
        newPath = path + "/" + name
        for heading in headings:
            value = "NA"
            for rHeading in reportHeadings:
                # Had some issues with the values being the wrong "type".
                # This if statement only accepts strings, floats, and integers
                if type(metadata[name][rHeading][heading]) is StringType or type(metadata[name][rHeading][heading]) \
                        is FloatType or type(metadata[name][rHeading][heading]) is IntType:
                    value = str(metadata[name][rHeading][heading])
            # Write the value to the spreadsheet
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
        jsonCollection = glob.glob("%s/reports/*Collection.json" % path)
        for jsonFile in jsonCollection:
            if os.path.isfile(jsonFile):
                os.remove(jsonFile)

        # Move assemblies and reports to appropriate Master repositories
        shutil.copyfile("%s/%s/%s_filteredAssembled.fasta" % (path, name, name),
                    "%s/AssemblyData/Assemblies/%s_filteredAssembled.fasta" % (repositoryPath, name))
        shutil.copyfile("%s/%s/%s_metadataReport.json" % (path, name, name),
                    "%s/AssemblyData/JsonReports/%s_metadataReport.json" % (repositoryPath, name))
    # Move the metadata spreadsheet
    combinedReport.close()
    shutil.copyfile("%s/%s_CombinedMetadataReport.tsv" % (reportPath, folderName),
                "%s/AssemblyData/SummaryReports/%s_CombinedMetadataReport.tsv" % (repositoryPath, folderName))


def functionsGoNOW(sampleNames, metadata, path):
    """Run the appropriate functions"""
    versionMetadata = versionIdentifyR.pipelineMetadata(path, metadata, sampleNames)
    reportWriter(sampleNames, versionMetadata, path)