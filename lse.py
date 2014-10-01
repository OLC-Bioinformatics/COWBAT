#!/usr/bin/env python
__author__ = 'akoziol'

from multiprocessing import Pool
# Shutil is useful for file moving/copying
import shutil
# Subprocess->call is used for making system calls
import subprocess
import sys
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# OS is used for file/folder manipulations
import os
import time
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob
import json

# referenceFile = glob.glob("*.fa*")
# references = ["%s/Best_Assemblies/%s" % (path, fastaFile) for fastaFile in referenceFile]
references = []
# inputData = {}

# targets = [reference.split('.')[0] for reference in referenceFile]
def referenceFiletoAssembly(path, sampleNames):
    for name in sampleNames:
        newPath = path + "/" + name
        references.append("%s/%s_filteredAssembled.fasta" % (newPath, name))


    # Create a dictionary of sorted tuples using zip
    inputData = dict(zip(references, sampleNames))
    return inputData
    # print json.dumps(inputData, sort_keys=True, indent=4)

def sampleFastq(path, sampleNames):
    """Sample the fastq files, so the processing doesn't take nearly as long"""
    for name in sampleNames:
        newPath = path + "/" + name
        if not os.path.isfile("%s/%s_R1_001_cor_sampled10000.fastq" % (newPath, name)):
            seqtkCall = "/home/blais/PycharmProjects/seqtk/seqtk sample -s seed=100 " \
                        "%s/%s_R1_001.cor.fastq 10000 > %s/%s_R1_001_cor_sampled10000.fastq " \
                        "&& /home/blais/PycharmProjects/seqtk/seqtk sample -s seed=100 " \
                        "%s/%s_R2_001.cor.fastq 10000 > %s/%s_R2_001_cor_sampled10000.fastq" \
                        % (newPath, name, newPath, name, newPath, name, newPath, name)
            # print seqtkCall
            os.system(seqtkCall)
            dotter()
        else:
            dotter()


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


dotcount = 0


def dotter():
    """This function is borrowed from Mike Knowles. It allows for the addition of pretty
    dots after every pool is finished its task. Additionally, it formats the dots such that
    there are only 80 dots per line, and the date is added at the start of each line"""
    global dotcount
    if dotcount <= 80:
        sys.stdout.write('.')
        # I added this flush command, as the dots were not being printed until the script
        # finished processing
        sys.stdout.flush()
        dotcount += 1
    else:
        sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
        dotcount = 0


def indexTargetsProcesses(path, inputData):
    sys.stdout.write('\nIndexing targets\n')
    indexTargetArgs = []
    scriptName = __name__
    if __name__ == scriptName:
        indexTargetsPool = Pool()
        # Initialise the pool of processes - it defaults to the number of processors
        for reference, target in inputData.iteritems():
            indexTargetArgs.append((reference, target, path))
        indexTargetsPool.map(indexTargets, indexTargetArgs)


def indexTargets((reference, target, path)):
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    filename = target.split('.')[0]
    newPath = path + "/" + target

    # Create a new path to be created (if necessary) for the generation of the range of k-mers
    indexPath = "%s/targets" % newPath
    # Call the make_path function to make folders as necessary
    make_path(indexPath)
    shutil.copy(reference, indexPath)
    # os.chdir(indexPath)
    indexFileSMI = "%s.smi" % filename
    if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
        indexCommand = "/bin/smalt index -k 20 -s 10 %s/%s_filteredAssembled %s" % (indexPath, target, reference)
        # print indexCommand
        subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # os.system(indexCommand)
        dotter()
    else:
        dotter()


def mappingProcesses(path, inputData):
    """Mapping threads!"""
    # os.chdir(path)
    print '\nPerforming reference mapping'
    mappingProcessesArgs = []
    scriptName = __name__
    if __name__ == scriptName:
        mappingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            mappingProcessesArgs.append((reference, target, path))
        mappingProcessesPool.map(mapping, mappingProcessesArgs)


def mapping((reference, target, path)):
    """Performs the mapping of the simulated reads to the targets"""
    filename = target.split('.')[0]
    newPath = path + "/" + target
    # os.chdir("%s/%s" % (path, target))
    fastq1 = "%s/%s_R1_001_cor_sampled10000.fastq" % (newPath, target)
    fastq2 = "%s/%s_R2_001_cor_sampled10000.fastq" % (newPath, target)
    # fastq1 = str("%s/%s" % (newPath, glob.glob("*R1_001_cor_sampled10000.fastq")[0]))
    # fastq2 = str("%s/%s" % (newPath, glob.glob("*R2_001_cor_sampled10000.fastq")[0]))
    filePath = "%s/tmp" % newPath
    # shutil.rmtree(filePath)
    make_path(filePath)
    targetPath = "%s/targets/%s" % (newPath, filename)
    if not os.path.isfile("%s/%s.bam" % (filePath, target)):
        smaltMap = "/bin/smalt map -o %s/%s.bam -f bam -n 24 -x %s_filteredAssembled %s %s" \
                   % (filePath, target, targetPath, fastq1, fastq2)
        # print smaltMap
        subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        # os.system(smaltMap)
        dotter()
    else:
        dotter()


def extractingProcesses(path, inputData):
    """Mapping threads!"""
    # os.chdir(path)
    print '\nExtracting insert sizes'
    extractingProcessesArgs = []
    scriptName = __name__
    if __name__ == scriptName:
        extractingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            extractingProcessesArgs.append((reference, target, path))
        extractingProcessesPool.map(extractInsertSize, extractingProcessesArgs)


def extractInsertSize((reference, target, path)):
    """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance between
    mapped paired reads"""
    # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
    newPath = path + "/" + target
    filePath = "%s/tmp" % newPath
    if not os.path.isfile("%s/%s_insertsizes.csv" % (filePath, target)):
        extractCommand = "samtools view %s/%s.bam | cut -f9 > %s/%s_insertsizes.csv" % (filePath, target, filePath, target)
        subprocess.call(extractCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        dotter()
    else:
        dotter()


def graphingProcesses(path, inputData):
    """Mapping threads!"""
    print '\nGraphing results'
    graphingProcessesArgs = []
    scriptName = __name__
    if __name__ == scriptName:
        graphingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            graphingProcessesArgs.append((target, path))
        graphingProcessesPool.map(graphing, graphingProcessesArgs)


def graphing((target, path)):
    """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance between
    mapped paired reads"""
    # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
    newPath = path + "/" + target
    filePath = "%s/tmp" % newPath
    newPath = "%s/insertSizes" % newPath
    make_path(newPath)
    os.chdir(newPath)
    if not os.path.isfile("%s/%s_insert_sizes.pdf" % (newPath, target)):
        graphingCommand = "Rscript /home/blais/PycharmProjects/LibrarySizeEstimator/insertsizes.R %s %s 1>/dev/null 2>/dev/null" % (filePath, target)
        # subprocess.call(graphingCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        os.system(graphingCommand)
        dotter()
    else:
        dotter()


def formatOutput(path, sampleNames, runTrimMetadata):
    insertSizePath = "%s/insertSizes" % path
    make_path(insertSizePath)
    # os.chdir("%s/insertSizes" % path)
    print("\nFormatting Outputs")
    # Determine the folder name by taking the last folder name from the path
    folderName = path.split('/')[-1]
    os.remove("%s/%s_insertSizes.csv" % (insertSizePath, folderName))
    # if not os.path.isfile("%s/%s_insertSizes.csv" % (insertSizePath, folderName)):
    with open("%s/%s_insertSizes.csv" % (insertSizePath, folderName), "a") as outputFile:
        outputFile.write("Strain\tMedian Insert Size\tStandard Deviation\n")
        for name in sampleNames:
            newPath = path + "/" + name
            if os.path.isfile("%s/insertSizes/%s_insert_sizes.txt" % (newPath, name)):
                infile = open("%s/insertSizes/%s_insert_sizes.txt" % (newPath, name), "r")
                inData = infile.read()
                data = inData.split("\t")
                infile.close()
                outputFile.write("%s\n" % inData)
                runTrimMetadata[name]["1.General"]["MedianInsertSize"] = data[1]
                runTrimMetadata[name]["1.General"]["InsertSizeStDev"] = data[2].rstrip()
                dotter()
            else:
                outputFile.write("%s\tN/A\tN/A\n" % name)
                runTrimMetadata[name]["1.General"]["MedianInsertSize"] = "N/A"
                runTrimMetadata[name]["1.General"]["InsertSizeStDev"] = "N/A"
                dotter()
    outputFile.close()
    return runTrimMetadata

def functionsGoNOW(sampleNames, path, runTrimMetadata):
    """Calls all the functions in a way that they can be multi-processed"""
    inputData = referenceFiletoAssembly(path, sampleNames)
    # print json.dumps(runTrimMetadata, sort_keys=True, indent=4)
    print "\nSampling fastq files."
    sampleFastq(path, sampleNames)
    indexTargetsProcesses(path, inputData)
    #Start the mapping operations
    mappingProcesses(path, inputData)
    extractingProcesses(path, inputData)
    graphingProcesses(path, inputData)
    os.chdir(path)
    runTrimInsertMetadata = formatOutput(path, sampleNames, runTrimMetadata)

    return runTrimInsertMetadata
    # print json.dumps(inputData, sort_keys=True, indent=4)

