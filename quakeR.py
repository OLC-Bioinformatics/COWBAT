__author__ = 'akoziol'

import os
from multiprocessing import Pool
import subprocess
import sys

# def quakePrep(sampleName, path):
#     for name in sampleName:
#         print name
#         os.chdir("%s/%s" % (path, name))


def quakePrepProcesses(sampleName, path):
    print 'Performing error correction on fastq files'
    quakePrepArgs = []
    if __name__ == 'quakeR':
        createQuakePool = Pool()
        # uses kmer, targets, readLength, foldCoverage
        for name in sampleName:
            quakePrepArgs.append((name, path))
        createQuakePool.map(runQuake, quakePrepArgs)
    # createVCFPool.terminate()
    # createVCFPool.join()


def runQuake((name, path)):
    """Creates the variant calling format files from which all relevant data can be pulled"""
    # print(name, path)
    forward = name + "_R1_001.fastq"
    reverse = name + "_R2_001.fastq"
    # print forward, reverse
    newPath = path + "/" + name
    # Check for the existence of the stats file - hopefully this will be created at the end
    # of the error correction processes
    if not os.path.isfile("%s/%s.stats.txt" % (newPath, forward)):
        # Create the list of fastq files required by quake
        fastqList = open("%s/%s_fastqFiles.txt" % (newPath, name), "wb")
        fastqList.write("%s/%s\t%s/%s" % (newPath, forward, newPath, reverse))
        fastqList.close()
        # Quake run command
        quakeRun = "/home/blais/Bioinformatics/Quake/bin/quake.py -f %s/%s_fastqFiles.txt -k 15 -p 1" % (newPath, name)
        # , stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb')
        subprocess.call(quakeRun, shell=True)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')






def functionsGoNOW(sampleNames, path):
    sampleNames = sampleNames
    path = path
    quakePrepProcesses(sampleNames, path)
