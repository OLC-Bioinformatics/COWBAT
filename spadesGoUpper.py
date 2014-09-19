__author__ = 'akoziol'

import os
from multiprocessing import Pool
import sys
import re
from Bio import SeqIO


def spadesPrepProcesses(sampleName, path):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    spadesPrepArgs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'spadesGoUpper':
        createSpadesPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            spadesPrepArgs.append((name, path))
        # This map function allows for multi-processing
        createSpadesPool.map(runSpades, spadesPrepArgs)


def runSpades((name, path)):
    """Performs necessary checks and runs SPAdes"""
    # Set up variables to keep commands clean looking
    scaffoldsFile = "scaffolds.fastg"
    newPath = path + "/" + name
    # Check for the existence of the scaffolds file - hopefully this will be created at the end of the run
    if not os.path.isfile("%s/spades_output/%s" % (newPath, scaffoldsFile)):
        # This is using a hard-coded path, as for some reason, when run within pycharm, spades.py could not
        # be located. Maybe the $PATH needs to be updated?
        # --continue
        forward = "%s/%s_R1_001.cor.fastq" %(newPath, name)
        reverse = "%s/%s_R2_001.cor.fastq" %(newPath, name)
        # There's an option to continue from checkpoints if the assembly is terminated prematurely within SPAdes,
        #  but the output directory must exist - if this directory exists, --continue, else don't --continue
        if os.path.isdir("%s/spades_output" % name):
            spadesRun = "/home/blais/Bioinformatics/SPAdes-3.1.1-Linux/bin/spades.py -k 21,33,55,77,99,127 " \
                        "--careful --continue --only-assembler --pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
        else:
            spadesRun = "/home/blais/Bioinformatics/SPAdes-3.1.1-Linux/bin/spades.py -k 21,33,55,77,99,127 --careful " \
                        "--only-assembler --pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
        # Run the command - subprocess.call would not run this command properly - no idea why - so using os.system instead
        # added 1>/dev/null to keep terminal output from being printed to screen
        os.system(spadesRun)
        # Print dots as per usual
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def contigFileFormatter(correctedFiles, path):
    """Changes the name of each contig from ">NODE_XXX_length..." to the name of the file ">OLC795_XXX..." """
    for name in correctedFiles:
        newPath = path + "/" + name
         #
        if os.path.isfile("%s/spades_output/contigs.fasta" % newPath) and not os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, name)):
            # http://biopython.org/wiki/SeqIO#Input.2FOutput_Example_-_Filtering_by_sequence_length
            over200bp = []
            for record in SeqIO.parse(open("%s/spades_output/contigs.fasta" % newPath, "rU"), "fasta"):
                if len(record.seq) >= 200:
                    # Add this record to our list
                    newID = re.sub("NODE", name, record.id)
                    record.id = newID
                    record.name = ''
                    record.description = ''
                    over200bp.append(record)
            # print "Found %s long sequences" % len(over200bp)
            formatted = open("%s/%s_filteredAssembled.fasta" % (newPath, name), "wb")
            SeqIO.write(over200bp, formatted, "fasta")
            formatted.close()


def functionsGoNOW(correctedFiles, path):
    """Run the helper function"""
    print("Assembling reads.")
    spadesPrepProcesses(correctedFiles, path)
    contigFileFormatter(correctedFiles, path)