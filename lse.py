#!/usr/bin/env python
# from multiprocessing import Pool
# Shutil is useful for file moving/copying
# import shutil
# Subprocess->call is used for making system calls
# import subprocess
# import sys
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
# import errno
# OS is used for file/folder manipulations
import os
# import time
# from collections import defaultdict
# import metadataFiller
import jsonReportR

__author__ = 'akoziol'
# # Initialise variables
# references = []
# dotcount = 0
#
#
# def make_dict():
#     """Makes Perl-style dictionaries"""
#     return defaultdict(make_dict)
#
#
# def referenceFiletoAssembly(path, sampleNames):
#     """Creates a dictionary of tuples of the reference genome and the assembled sequences"""
#     for name in sampleNames:
#         newPath = path + "/" + name
#         references.append("%s/%s_filteredAssembled.fasta" % (newPath, name))
#     # Create a dictionary of sorted tuples using zip
#     inputData = dict(zip(references, sampleNames))
#     return inputData
#
#
# def sampleFastq(path, sampleNames, metadata, commands):
#     """Sample the fastq files, so the processing doesn't take nearly as long"""
#     for name in sampleNames:
#         newPath = path + "/" + name
#         # Randomly samples 10 000 reads with a seed of 100
#         # if not os.path.isfile("%s/%s_R1_001_sampled10000.fastq" % (newPath, name)) and not \
#         #     os.path.isfile("%s/%s_R1_001_sampled10000.fastq.xz" % (newPath, name)):
#         if not commands[name]["seqtkSamplingCommand"]:
#             seqtkCall = "seqtk sample -s seed=100 %s/%s_R1_001.cor.fastq 10000 > %s/%s_R1_001_sampled10000.fastq " \
#                         "&& seqtk sample -s seed=100 %s/%s_R2_001.cor.fastq 10000 > %s/%s_R2_001_sampled10000.fastq" \
#                         % (newPath, name, newPath, name, newPath, name, newPath, name)
#             os.system(seqtkCall)
#             metadata[name]["7.PipelineCommands"]["seqtkSamplingCommand"] = seqtkCall
#             dotter()
#         else:
#             dotter()
#             metadata[name]["7.PipelineCommands"]["seqtkSamplingCommand"] = commands[name]["seqtkSamplingCommand"]
#     return metadata
#
#
# def make_path(inPath):
#     """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
#     does what is indicated by the URL"""
#     try:
#         os.makedirs(inPath)
#         os.chmod(inPath, 0777)
#     except OSError as exception:
#         if exception.errno != errno.EEXIST:
#             raise
#
#
# def dotter():
#     """This function is borrowed from Mike Knowles. It allows for the addition of pretty
#     dots after every pool is finished its task. Additionally, it formats the dots such that
#     there are only 80 dots per line, and the date is added at the start of each line"""
#     global dotcount
#     if dotcount <= 80:
#         sys.stdout.write('.')
#         # I added this flush command, as the dots were not being printed until the script
#         # finished processing
#         sys.stdout.flush()
#         dotcount += 1
#     else:
#         sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
#         dotcount = 0
#
#
# def indexTargetsProcesses(path, inputData, metadata, commands):
#     """Allows for multiprocessing of smalt index on targets"""
#     sys.stdout.write('\nIndexing targets\n')
#     indexTargetArgs = []
#     output = {}
#     scriptName = __name__
#     if __name__ == scriptName:
#         indexTargetsPool = Pool()
#         # Initialise the pool of processes - it defaults to the number of processors
#         for reference, target in inputData.iteritems():
#             indexTargetArgs.append((reference, target, path, metadata, commands))
#         output = indexTargetsPool.map(indexTargets, indexTargetArgs)
#     return output
#
#
# def indexTargets((reference, target, path, metadata, commands)):
#     """Performs smalt index on the targets"""
#     newPath = path + "/" + target
#     # if os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, target)):
#     if not commands[target]["smaltIndexCommand"]:
#         filename = target.split('.')[0]
#         # Create a new path to be created (if necessary) for the generation of the range of k-mers
#         indexPath = "%s/targets" % newPath
#         # Call the make_path function to make folders as necessary
#         make_path(indexPath)
#         shutil.copy(reference, indexPath)
#         indexFileSMI = "%s.smi" % filename
#         # if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
#         indexCommand = "smalt index -k 20 -s 10 %s/%s_filteredAssembled %s" % (indexPath, target, reference)
#         subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
#         metadata[target]["7.PipelineCommands"]["smaltIndexCommand"] = indexCommand
#         dotter()
#         # else:
#         #     dotter()
#         #     metadata[target]["7.Pipeline"]["smaltIndexCommand"] = "N/A"
#     else:
#         dotter()
#         metadata[target]["7.PipelineCommands"]["smaltIndexCommand"] = commands[target]["smaltIndexCommand"]
#     return metadata
#
#
# def mappingProcesses(path, inputData, metadata, commands):
#     """Mapping threads!"""
#     print '\nPerforming reference mapping'
#     mappingProcessesArgs = []
#     output = {}
#     scriptName = __name__
#     if __name__ == scriptName:
#         mappingProcessesPool = Pool()
#         # uses target
#         for reference, target in inputData.iteritems():
#             mappingProcessesArgs.append((target, path, metadata, commands))
#         output = mappingProcessesPool.map(mapping, mappingProcessesArgs)
#     return output
#
#
# def mapping((target, path, metadata, commands)):
#     """Performs the mapping of the sampled reads to the targets"""
#     newPath = path + "/" + target
#     # if os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, target)):
#     if not commands[target]["smaltMapCommand"]:
#         filename = target.split('.')[0]
#         fastq1 = "%s/%s_R1_001_sampled10000.fastq" % (newPath, target)
#         fastq2 = "%s/%s_R2_001_sampled10000.fastq" % (newPath, target)
#         filePath = "%s/tmp" % newPath
#         make_path(filePath)
#         targetPath = "%s/targets/%s" % (newPath, filename)
#         # if not os.path.isfile("%s/%s.bam" % (filePath, target)):
#         smaltMap = "smalt map -o %s/%s.bam -f bam -n 24 -x %s_filteredAssembled %s %s" \
#                    % (filePath, target, targetPath, fastq1, fastq2)
#         subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
#         metadata[target]["7.PipelineCommands"]["smaltMapCommand"] = smaltMap
#         dotter()
#         # else:
#         #     metadata[target]["7.Pipeline"]["SmaltMapCommand"] = "N/A"
#         #     dotter()
#     else:
#         metadata[target]["7.PipelineCommands"]["smaltMapCommand"] = commands[target]["smaltMapCommand"]
#         dotter()
#     return metadata
#
#
# def extractingProcesses(path, inputData, metadata, commands):
#     """Mapping threads!"""
#     print '\nExtracting insert sizes'
#     extractingProcessesArgs = []
#     output = {}
#     scriptName = __name__
#     if __name__ == scriptName:
#         extractingProcessesPool = Pool()
#         # uses target
#         for reference, target in inputData.iteritems():
#             extractingProcessesArgs.append((target, path, metadata, commands))
#         output = extractingProcessesPool.map(extractInsertSize, extractingProcessesArgs)
#     return output
#
# def extractInsertSize((target, path, metadata, commands)):
#     """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance
# between
#     mapped paired reads"""
#     # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
#     fileSize = 0
#     newPath = path + "/" + target
#     filePath = "%s/tmp" % newPath
#     # if os.path.isfile("%s/%s_insertsizes.csv" % (filePath, target)):
#     #     fileSize = os.stat("%s/%s_insertsizes.csv" % (filePath, target)).st_size
#     if not commands[target]["samtoolsCommand"]:
#     # if os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, target)) and fileSize == 0:
#         extractCommand = "samtools view %s/%s.bam | cut -f9 > %s/%s_insertsizes.csv" % (filePath, target, filePath,
# target)
#         subprocess.call(extractCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
#         metadata[target]["7.PipelineCommands"]["samtoolsCommand"] = extractCommand
#         dotter()
#     else:
#         metadata[target]["7.PipelineCommands"]["samtoolsCommand"] = commands[target]["samtoolsCommand"]
#         dotter()
#     return metadata
#
#
# def graphingProcesses(path, inputData, metadata, commands):
#     """Mapping threads!"""
#     print '\nGraphing results'
#     graphingProcessesArgs = []
#     output = {}
#     scriptName = __name__
#     if __name__ == scriptName:
#         graphingProcessesPool = Pool()
#         # uses target
#         for reference, target in inputData.iteritems():
#             graphingProcessesArgs.append((target, path, metadata, commands))
#         output = graphingProcessesPool.map(graphing, graphingProcessesArgs)
#     return output
#
#
# def graphing((target, path, metadata, commands)):
#     """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance
# between
#     mapped paired reads"""
#     newPath = path + "/" + target
#     # if os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, target)):
#     if not commands[target]["InsertSizeGraphingCommand"]:
#         filePath = "%s/tmp" % newPath
#         newPath = "%s/insertSizes" % newPath
#         make_path(newPath)
#         os.chdir(newPath)
#         # if not os.path.isfile("%s/%s_insert_sizes.pdf" % (newPath, target)):
#             #  1>/dev/null 2>/dev/null
#         graphingCommand = "insertsizes.R %s %s 1>/dev/null 2>/dev/null" % (filePath, target)
#         os.system(graphingCommand)
#         metadata[target]["7.PipelineCommands"]["InsertSizeGraphingCommand"] = graphingCommand
#         dotter()
#         # else:
#         #     dotter()
#         #     metadata[target]["7.Pipeline"]["InsertSizeGraphingCommand"] = "N/A"
#     else:
#         dotter()
#         metadata[target]["7.PipelineCommands"]["InsertSizeGraphingCommand"] =
# commands[target]["InsertSizeGraphingCommand"]
#     return metadata
#
#
# def formatOutput(path, sampleNames, runTrimMetadata):
#     """Gets all the insert size metrics into a single report, as well as adding those
#     data to the metadata dictionary"""
#     insertSizePath = "%s/insertSizes" % path
#     make_path(insertSizePath)
#     print("\nFormatting Outputs")
#     # Determine the folder name by taking the last folder name from the path
#     folderName = path.split('/')[-1]
#     # As the file is opened to append - it must be deleted each time through the pipeline
#     if os.path.isfile("%s/%s_insertSizes.csv" % (insertSizePath, folderName)):
#         os.remove("%s/%s_insertSizes.csv" % (insertSizePath, folderName))
#     with open("%s/%s_insertSizes.csv" % (insertSizePath, folderName), "a") as outputFile:
#         outputFile.write("Strain\tMean Insert Size\tStandard Deviation\n")
#         for name in sampleNames:
#             newPath = path + "/" + name
#             if os.path.isfile("%s/insertSizes/%s_insert_sizes.txt" % (newPath, name)):
#                 infile = open("%s/insertSizes/%s_insert_sizes.txt" % (newPath, name), "r")
#                 inData = infile.read()
#                 data = inData.split("\t")
#                 infile.close()
#                 outputFile.write("%s\n" % inData)
#                 insertSize = "%.2f" % float(data[1])
#                 runTrimMetadata[name]["1.General"]["MeanInsertSize"] = insertSize
#                 runTrimMetadata[name]["1.General"]["InsertSizeStDev"] = float(data[2].rstrip())
#                 dotter()
#             else:
#                 outputFile.write("%s\tN/A\tN/A\n" % name)
#                 runTrimMetadata[name]["1.General"]["MeanInsertSize"] = "N/A"
#                 runTrimMetadata[name]["1.General"]["InsertSizeStDev"] = "N/A"
#                 dotter()
#     outputFile.close()
#     return runTrimMetadata


def insertsize(path, samplenames, metadata, numreads):
    """
    Extracts the insert size and its deviation from the spades.log file
    :param numreads: the number of reads (1 for single reads, 2 for paired-end)
    :param metadata: dictionary of metadata
    :param samplenames: a list of sample names
    :param path: the path of the analyses
    """
    for name in samplenames:
        newpath = path + "/" + name
        # Set the name of the log file
        spadeslogfile = '{}/spades_output/spades.log'.format(newpath)
        # Only look if the spades output folder exists, and if there are two fastq files (can't find the insert
        # size of single reads
        if os.path.isfile(spadeslogfile) and numreads == 2:
            # Open the log file
            with open(spadeslogfile, 'rb') as spadeslog:
                # Iterate through the file
                for line in spadeslog:
                    # Find the line with the insert size on it. Will look something like this:
                    """
                    0:02:07.605   144M / 9G    INFO    General (pair_info_count.cpp : 191) \
                    Insert size = 240.514, deviation = 105.257, left quantile = 142, right quantile = 384, \
                    read length = 301
                    """
                    if 'Insert size =' in line:
                        # Extract the relevant data and add it to the metadata
                        metadata[name]["1.General"]["MeanInsertSize"] = line.split('= ')[1].split(',')[0]
                        metadata[name]["1.General"]["InsertSizeStDev"] = line.split('= ')[2].split(',')[0]
        # Otherwise populate with NA
        else:
            metadata[name]["1.General"]["MeanInsertSize"] = 'NA'
            metadata[name]["1.General"]["InsertSizeStDev"] = 'NA'
    return metadata


def functionsgonow(samplenames, path, metadata, numreads):
    runtriminsertmetadata = insertsize(path, samplenames, metadata, numreads)
#     """Calls all the functions in a way that they can be multi-processed"""
#     inputData = referenceFiletoAssembly(path, sampleNames)
#     print "\nSampling fastq files."
#     sampleMeta = sampleFastq(path, sampleNames, runTrimMetadata, commands)
#     indexList = indexTargetsProcesses(path, inputData, sampleMeta, commands)
#     indexMeta = metadataFiller.filler(runTrimMetadata, indexList)
#     #Start the mapping operations
#     mappingList = mappingProcesses(path, inputData, indexMeta, commands)
#     mappingMeta = metadataFiller.filler(runTrimMetadata, mappingList)
#     extractingList = extractingProcesses(path, inputData, mappingMeta, commands)
#     extractingMeta = metadataFiller.filler(runTrimMetadata, extractingList)
#     graphingList = graphingProcesses(path, inputData, extractingMeta, commands)
#     graphingMeta = metadataFiller.filler(runTrimMetadata, graphingList)
#     os.chdir(path)
#     runTrimInsertMetadata = formatOutput(path, sampleNames, graphingMeta)
    jsonReportR.jsonR(samplenames, path, runtriminsertmetadata, "Collection")
    return runtriminsertmetadata
