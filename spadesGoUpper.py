__author__ = 'akoziol'

import os
from multiprocessing import Pool
import sys
import re
from Bio import SeqIO
import shutil
import errno
import subprocess
import metadataFiller
# add spades.py to PYTHONPATH and interact directly
# import spades
import json
from collections import defaultdict
import glob
import jsonReportR

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

flagMetadata = defaultdict(make_dict)


def spadesPrepProcesses(sampleName, path, fLength, metadata, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    spadesPrepArgs = []
    flagMetadata = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'spadesGoUpper':
        createSpadesPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            spadesPrepArgs.append((name, path, fLength, metadata, commands))
        # This map function allows for multi-processing
        flagMetadata = createSpadesPool.map(runSpades, spadesPrepArgs)
    #flagMetadata = metadata
    # for item in flagMetadata:
    return flagMetadata


def runSpades((name, path, fLength, metadata, commands)):
    """Performs necessary checks and runs SPAdes"""
    # Set up variables to keep commands clean looking
    # contigsFile = "contigs.fasta"
    newPath = path + "/" + name
    forward = ""
    reverse = ""
    # Check for the existence of the scaffolds file - hopefully this will be created at the end of the run
    if not os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, name)) or not commands[name]["SPAdesCommand"]:
        # # As this now re-runs previously assembled strains lacking the appropriate command in the metadata, the output file must
        # # be removed prior to
        # if os.path.isdir("%s/%s/spades_output" % (path, name)):
        #     shutil.rmtree("%s/%s/spades_output" % (path, name))
        forward = "%s/%s_R1_001.cor.fastq" % (newPath, name)
        reverse = "%s/%s_R2_001.cor.fastq" % (newPath, name)
        # forwardFile = glob.glob("%s/*1.cor.fastq" % newPath)
        # if forwardFile:
        #     forward = forwardFile[0]
        # reverseFile = glob.glob("%s/*2.cor.fastq" % newPath)
        # if reverseFile:
        #     reverse = reverseFile[0]
        # forward = "%s/%s_R1_001.fastq" % (newPath, name)
        # reverse = "%s/%s_R2_001.fastq" % (newPath, name)
        # There's an option to continue from checkpoints if the assembly is terminated prematurely within SPAdes,
        #  but the output directory must exist - if this directory exists, --continue, else don't --continue
        # /home/blais/Bioinformatics/SPAdes-3.1.1-Linux/bin/
        if fLength > 50:
            if os.path.isdir("%s/spades_output" % name):
                spadesRun = "spades.py -k 21,33,55,77,99,127 --careful --continue " \
                            "--only-assembler --pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
                metadata[name]["7.PipelineCommands"]["SPAdesCommand"] = spadesRun
            else:
                spadesRun = "spades.py -k 21,33,55,77,99,127 --careful --only-assembler " \
                            "--pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
                metadata[name]["7.PipelineCommands"]["SPAdesCommand"] = spadesRun
        else:
            if os.path.isdir("%s/spades_output" % name):
                spadesRun = "spades.py -k 21,33,55,77,99,127 --careful --continue " \
                            "--only-assembler --s1 %s -o %s/spades_output 1>/dev/null" % (reverse, newPath)
                metadata[name]["7.PipelineCommands"]["SPAdesCommand"] = spadesRun
            else:
                spadesRun = "spades.py -k 21,33,55,77,99,127 --careful --only-assembler " \
                            "--s1 %s -o %s/spades_output 1>/dev/null" % (reverse, newPath)
                metadata[name]["7.PipelineCommands"]["SPAdesCommand"] = spadesRun
        # Run the command - subprocess.call would not run this command properly - no idea why - so using os.system instead
        # added 1>/dev/null to keep terminal output from being printed to screen
        os.system(spadesRun)
        # Print dots as per usual
        sys.stdout.write('.')
        return metadata
    else:
        sys.stdout.write('.')
        metadata[name]["7.PipelineCommands"]["SPAdesCommand"] = commands[name]["SPAdesCommand"]
        return metadata


def contigFileFormatter(correctedFiles, path, metadata):
    """Changes the name of each contig from ">NODE_XXX_length..." to the name of the file ">OLC795_XXX..." """
    for name in correctedFiles:
        newPath = path + "/" + name
        # Ensures that the contigs file is present, but the renamed, manipulated file is not
        # os.remove("%s/%s_filteredAssembled.fasta" % (newPath, name))
        # and not os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, name))
        if os.path.isfile("%s/spades_output/contigs.fasta" % newPath):
            # http://biopython.org/wiki/SeqIO#Input.2FOutput_Example_-_Filtering_by_sequence_length
            over1000bp = []
            lengthCov = 0
            for record in SeqIO.parse(open("%s/spades_output/contigs.fasta" % newPath, "rU"), "fasta"):
                cov = float(record.id.split("_")[5])
                # Include only contigs greater than 1000 bp in length
                 #and cov > 10
                if len(record.seq) >= 1000:
                    # Add this record to our list
                    newID = re.sub("NODE", name, record.id)
                    lengthCov += (float(newID.split("_")[-5]) * float(newID.split("_")[-3]))
                    record.id = newID
                    record.name = ''
                    record.description = ''
                    over1000bp.append(record)
            metadata[name]["2.Assembly"]["totalBasesxCoverage"] = lengthCov
            # print metadata[name]["2.Assembly"]["totalBasesxCoverage"]
            # print "Found %s long sequences" % len(over200bp)
            fileName = "%s/%s_filteredAssembled.fasta" % (newPath, name)
            formatted = open(fileName, "wb")
            SeqIO.write(over1000bp, formatted, "fasta")
            formatted.close()
            # Move the files to a BestAssemblies folder
            assemblyPath = "%s/BestAssemblies" % path
            make_path(assemblyPath)
            if not os.path.isfile("%s/%s" % (assemblyPath, fileName)):
                shutil.copy(fileName, assemblyPath)
        fileName = "%s/%s_filteredAssembled.fasta" % (newPath, name)
        lengthCov = 0
        if os.path.isfile(fileName):
            for record in SeqIO.parse(open(fileName, "rU"), "fasta"):
                # lengthCov = 0
                newID = re.sub("NODE", name, record.id)
                lengthCov += (float(newID.split("_")[-5]) * float(newID.split("_")[-3]))
            metadata[name]["2.Assembly"]["totalBasesxCoverage"] = lengthCov
    return metadata


def completionist(correctedFiles, path):
    """Function to determine for which strains assembly was successful. Creates a list of the successful
     assemblies that will be further processed by the pipeline"""
    # Initialise the list
    assembledFiles = []
    for name in correctedFiles:
        # Check if the assembly has completed successfully
        if os.path.isfile("%s/%s/%s_filteredAssembled.fasta" % (path, name, name)):
            # Append successful assemblies to the list
            assembledFiles.append(name)
        # Else inform that assembly failed
        else:
            print "%s could not be assembled!" % name
    return assembledFiles


def pipelineMetadata(path, metadata, sampleNames):
    # Update these values as required
    #Quast
    quastOutput = subprocess.Popen(["quast.py"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = quastOutput.communicate()
    quastVer = err.split("\n")[1].strip("Version ").replace(" ,", ",").replace(",", ";")

    # SMALT
    smaltOutput = subprocess.Popen(["smalt"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = smaltOutput.communicate()
    smaltVer = out.split("\n")[2].strip(" ").strip("(").strip(")").strip("version: ")

    # Samtools
    samOutput = subprocess.Popen(["samtools"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = samOutput.communicate()
    samVer = err.split("\n")[2].strip("Version: ")

    # Blast
    blastnOutput = subprocess.Popen(["blastn -version"], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
    out, err = blastnOutput.communicate()
    blastVer = out.split("\n")[0].strip("blastn: ")

    for name in sampleNames:

        with open("%s/%s/spades_output/spades.log" % (path, name)) as spadesLog:
            for line in spadesLog:
                if "SPAdes version" in line:
                    spadesVer = line.strip("SPAdes version: ").rstrip()
                elif "Python version" in line:
                    pythonVer = line.strip("Python version: ").rstrip()
                elif "OS" in line:
                    OS = line.strip("OS: ").rstrip()

        quakeVer = "0.3"
        commit = "b737e2c52f59c541062a5f71c5e08087ee238c1e 2015-02-17 Adam Koziol"

        if not os.path.isfile("%s/%s/%s_programVersions.tsv" % (path, name, name)):

            versions = open("%s/%s/%s_programVersions.tsv" % (path, name, name), "wb")
            versions.write("spadesVersion\tquastVersion\tquakeVersion\tSmaltVersion\tSamtools\tBlastVersion\t"
                       "PythonVersion\tOS\tCommit\n")
            versions.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (spadesVer, quastVer, quakeVer, smaltVer, samVer, blastVer, pythonVer, OS, commit))

            # print name, quastVer, smaltVer, samVer, blastVer, spadesVer, pythonVer, OS
            metadata[name]["8.PipelineVersions"]["SPAdesVersion"] = spadesVer
            metadata[name]["8.PipelineVersions"]["QUASTVersion"] = quastVer
            metadata[name]["8.PipelineVersions"]["QuakeVersion"] = quakeVer
            metadata[name]["8.PipelineVersions"]["SmaltVersion"] = smaltVer
            metadata[name]["8.PipelineVersions"]["SamtoolsVersion"] = samVer
            metadata[name]["8.PipelineVersions"]["BlastVersion"] = blastVer
            metadata[name]["8.PipelineVersions"]["PythonVersion"] = pythonVer
            metadata[name]["8.PipelineVersions"]["OS"] = OS
            metadata[name]["8.PipelineVersions"]["PipelineVersion"] = commit
        else:
            pipelineVersions = open("%s/%s/%s_programVersions.tsv" % (path, name, name)).readlines()[1]
            versions = pipelineVersions.split("\t")
            metadata[name]["8.PipelineVersions"]["SPAdesVersion"] = versions[0]
            metadata[name]["8.PipelineVersions"]["QUASTVersion"] = versions[1]
            metadata[name]["8.PipelineVersions"]["QuakeVersion"] = versions[2]
            metadata[name]["8.PipelineVersions"]["SmaltVersion"] = versions[3]
            metadata[name]["8.PipelineVersions"]["SamtoolsVersion"] = versions[4]
            metadata[name]["8.PipelineVersions"]["BlastVersion"] = versions[5]
            metadata[name]["8.PipelineVersions"]["PythonVersion"] = versions[6]
            metadata[name]["8.PipelineVersions"]["OS"] = versions[7]
            metadata[name]["8.PipelineVersions"]["PipelineVersion"] = versions[8].rstrip()
    return metadata


def functionsGoNOW(correctedFiles, path, metadata, fLength, commands):
    """Run the helper function"""
    print("\nAssembling reads.")
    flagMetadataList = spadesPrepProcesses(correctedFiles, path, fLength, metadata, commands)
    flagMetadata = metadataFiller.filler(metadata, flagMetadataList)
    updatedMetadata = contigFileFormatter(correctedFiles, path, flagMetadata)
    assembledFiles = completionist(correctedFiles, path)
    moreMetadata = pipelineMetadata(path, updatedMetadata, assembledFiles)
    # print json.dumps(moreMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    jsonReportR.jsonR(correctedFiles, path, moreMetadata, "Collection")
    return moreMetadata, assembledFiles