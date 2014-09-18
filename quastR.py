__author__ = 'akoziol'

import os
import glob
from multiprocessing import Pool
import sys


def performQuast(newPath, quastCall):
    """Runs quast"""
    # Check to see if the file "report.txt" exists in the quast_results directory
    # this hopefully indicates that quast ran properly
    if os.path.isfile("%s/quast_results/report.tex" % newPath):
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
    # Check to see if there is a directory named referenceGenome - if there is, use the reference genome in that folder
    # in the quast analyses. If not, then perform quast without a reference genome.
    if os.path.isdir("%s/referenceGenome" % newPath):
        referenceGenome = glob.glob("%s/referenceGenome/*" % newPath)
        quastCall = "/home/blais/Bioinformatics/quast-2.3/quast.py -R %s --threads 24 --gage --scaffolds %s/spades_output/scaffolds.fasta -o %s/quast_results 1>/dev/null" % (referenceGenome[0], newPath, newPath)
        performQuast(newPath, quastCall)
    else:
        quastCall = "/home/blais/Bioinformatics/quast-2.3/quast.py --threads 24 --scaffolds %s/spades_output/scaffolds.fasta -o %s/quast_results 1>/dev/null" % (newPath, newPath)
        performQuast(newPath, quastCall)


def functionsGoNOW(sampleNames, path):
    print "\nPerforming quality checks on assemblies."
    quastProcesses(sampleNames, path)
