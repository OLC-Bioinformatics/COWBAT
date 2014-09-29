__author__ = 'akoziol'

import os
import glob
from multiprocessing import Pool
import sys
import re
import time
import errno
import json
import shutil


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def performQuast(newPath, quastCall):
    """Runs quast"""
    # Check to see if the file "report.txt" exists in the quast_results directory
    # this hopefully indicates that quast ran properly
    if os.path.isfile("%s/quast_results/report.tsv" % newPath):
        sys.stdout.write('.')
    else:
        # Run the command
        os.system(quastCall)
        sys.stdout.write('.')


def quastProcesses(sampleNames, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to quast
    assembly metrics processing"""
    quastPrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'quastR':
        createQuastPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleNames:
            quastPrepArgs.append((name, path))
        # This map function allows for multi-processing
        createQuastPool.map(quasting, quastPrepArgs)


def quasting((name, path)):
    """Performs quast analysis on the assemblies"""
    newPath = path + "/" + name
    # make_path("%s/quast_results" % newPath)
    # shutil.rmtree("%s/quast_results" % newPath)
    # Check to see if there is a directory named referenceGenome - if there is, use the reference genome in that folder
    # in the quast analyses. If not, then perform quast without a reference genome.
    if os.path.isdir("%s/referenceGenome" % newPath):
        referenceGenome = glob.glob("%s/referenceGenome/*" % newPath)
        quastCall = "/home/blais/Bioinformatics/quast-2.3/quast.py -R %s --threads 24 --gage %s/%s_filteredAssembled.fasta -o %s/quast_results 1>/dev/null" % (referenceGenome[0], newPath, name, newPath)
        performQuast(newPath, quastCall)
    else:
        quastCall = "/home/blais/Bioinformatics/quast-2.3/quast.py --threads 24 %s/spades_output/contigs.fasta -o %s/quast_results 1>/dev/null" % (newPath, newPath)
        performQuast(newPath, quastCall)


def quastMetadata(sampleNames, path, runTrimMetadata):
    """Pulls the metadata from the quast assembly reports"""
    for name in sampleNames:
        newPath = path + "/" + name
        if os.path.isdir("%s/referenceGenome" % newPath):
            referenceGenome = glob.glob("%s/referenceGenome/*" % newPath)
        # Populate the dictionary with the referenceGenome - this may be moved to the rMLST module
            runTrimMetadata[name]["1.General"]["referenceGenome"] = re.split("\.", re.split("/", referenceGenome[0])[-1])[0]
        # This is just here as a placeholder until the rMLST module is functional
        runTrimMetadata[name]["1.General"]["rMLST_sequenceType"] = "N/A"
        runTrimMetadata[name]["1.General"]["fileName"] = name
        runTrimMetadata[name]["1.General"]["filelocation"] = newPath
        runTrimMetadata[name]["1.General"]["assemblyDate"]= time.strftime("%Y-%m-%d")
        runTrimMetadata[name]["2.Assembly"]["kmerRange"] = "21,33,55,77,99,127"
        runTrimMetadata[name]["2.Assembly"]["assemblyType"] = "PE"
        if not os.path.isfile("%s/quast_results/report.tsv" % newPath):
            print "There was an issue getting the metadata from %s" % name
        else:
            report = open("%s/quast_results/transposed_report.tsv" % newPath, "r")
            for line in report:
                if re.search("Assembly", line):
                    pass
                else:
                    if os.path.isfile("%s/quast_results/gage_report.tsv" % newPath):
                        subline = line.split("\t")
                        # Populate the dictionary
                        runTrimMetadata[name]["2.Assembly"]["Assembly"] = subline[0]
                        runTrimMetadata[name]["2.Assembly"]["NumContigs"] = subline[1]
                        runTrimMetadata[name]["2.Assembly"]["NumContigsOver1000bp"] = subline[2]
                        runTrimMetadata[name]["2.Assembly"]["TotalLength"] = subline[3]
                        runTrimMetadata[name]["2.Assembly"]["TotalLengthOver1000bp"] = subline[4]
                        runTrimMetadata[name]["2.Assembly"]["NumContigsOver500bp"] = subline[5]
                        runTrimMetadata[name]["2.Assembly"]["LargestContig"] = subline[6]
                        runTrimMetadata[name]["2.Assembly"]["TotalLengthOver500bp"] = subline[7]
                        runTrimMetadata[name]["2.Assembly"]["ReferenceLength"] = subline[8]
                        runTrimMetadata[name]["2.Assembly"]["percentGC"] = subline[9]
                        runTrimMetadata[name]["2.Assembly"]["ReferencePercentGC"] = subline[10]
                        runTrimMetadata[name]["2.Assembly"]["N50"] = subline[11]
                        runTrimMetadata[name]["2.Assembly"]["NG50"] = subline[12]
                        runTrimMetadata[name]["2.Assembly"]["N75"] = subline[13]
                        runTrimMetadata[name]["2.Assembly"]["NG75"] = subline[14]
                        runTrimMetadata[name]["2.Assembly"]["L50"] = subline[15]
                        runTrimMetadata[name]["2.Assembly"]["LG50"] = subline[16]
                        runTrimMetadata[name]["2.Assembly"]["L75"] = subline[17]
                        runTrimMetadata[name]["2.Assembly"]["LG75"] = subline[18]
                        runTrimMetadata[name]["2.Assembly"]["NumMisassemblies"] = subline[19]
                        runTrimMetadata[name]["2.Assembly"]["NumMisassembledContigs"] = subline[20]
                        runTrimMetadata[name]["2.Assembly"]["MisassembledContigsLength"] = subline[21]
                        runTrimMetadata[name]["2.Assembly"]["NumLocalMisassemblies"] = subline[22]
                        runTrimMetadata[name]["2.Assembly"]["NumUnalignedContigs"] = subline[23]
                        runTrimMetadata[name]["2.Assembly"]["UnalignedLength"] = subline[24]
                        runTrimMetadata[name]["2.Assembly"]["percentGenomeFraction"] = subline[25]
                        runTrimMetadata[name]["2.Assembly"]["DuplicationRatio"] = subline[26]
                        runTrimMetadata[name]["2.Assembly"]["NumNsPer100kbp"] = subline[27]
                        runTrimMetadata[name]["2.Assembly"]["NumMismatchesPer100kbp"] = subline[28]
                        runTrimMetadata[name]["2.Assembly"]["NumIndelsPer100kbp"] = subline[29]
                        runTrimMetadata[name]["2.Assembly"]["LargestAlignment"] = subline[30].rstrip()
                    else:
                        # As above, but since the gage analysis wasn't performed,
                        # populate the dictionary with N/A where appropriate
                        subline = line.split("\t")
                        runTrimMetadata[name]["2.Assembly"]["Assembly"] = subline[0]
                        runTrimMetadata[name]["2.Assembly"]["NumContigs"] = subline[1]
                        runTrimMetadata[name]["2.Assembly"]["NumContigsOver1000bp"] = subline[2]
                        runTrimMetadata[name]["2.Assembly"]["TotalLength"] = subline[3]
                        runTrimMetadata[name]["2.Assembly"]["TotalLengthOver1000bp"] = subline[4]
                        runTrimMetadata[name]["2.Assembly"]["NumContigsOver500bp"] = subline[5]
                        runTrimMetadata[name]["2.Assembly"]["LargestContig"] = subline[6]
                        runTrimMetadata[name]["2.Assembly"]["TotalLengthOver500bp"] = subline[7]
                        runTrimMetadata[name]["2.Assembly"]["ReferenceLength"] = subline[8]
                        runTrimMetadata[name]["2.Assembly"]["percentGC"] = subline[9]
                        runTrimMetadata[name]["2.Assembly"]["N50"] = subline[10]
                        runTrimMetadata[name]["2.Assembly"]["N75"] = subline[11]
                        runTrimMetadata[name]["2.Assembly"]["L50"] = subline[12]
                        runTrimMetadata[name]["2.Assembly"]["L75"] = subline[13]
                        # These values aren't determined in a non-gage analysis, so in order to keep
                        # the schema consistent they are populated with "N/A"
                        runTrimMetadata[name]["2.Assembly"]["ReferencePercentGC"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NG50"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NG75"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["LG50"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["LG75"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumMisassemblies"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumMisassembledContigs"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["MisassembledContigsLength"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumLocalMisassemblies"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumUnalignedContigs"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["UnalignedLength"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["percentGenomeFraction"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["DuplicationRatio"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumNsPer100kbp"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumMismatchesPer100kbp"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["NumIndelsPer100kbp"] = "N/A"
                        runTrimMetadata[name]["2.Assembly"]["LargestAlignment"] = "N/A"
    return runTrimMetadata


def functionsGoNOW(sampleNames, path, runTrimMetadata):
    print "\nPerforming quality checks on assemblies."
    quastProcesses(sampleNames, path)
    print "\nCollating metadata."
    runTrimAssemblyMetadata = quastMetadata(sampleNames, path, runTrimMetadata)
    return runTrimAssemblyMetadata
    # print json.dumps(runTrimAssemblyMetadata, sort_keys=True, indent=4)
