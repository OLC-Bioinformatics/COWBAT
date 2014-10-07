__author__ = 'akoziol'

import os
from multiprocessing import Pool
import sys
import re
from Bio import SeqIO
import shutil
import errno


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


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
    contigsFile = "contigs.fasta"
    newPath = path + "/" + name
    # Check for the existence of the scaffolds file - hopefully this will be created at the end of the run
    if not os.path.isfile("%s/spades_output/%s" % (newPath, contigsFile)):
        # This is using a hard-coded path, as for some reason, when run within pycharm, spades.py could not
        # be located. Maybe the $PATH needs to be updated?
        # --continue
        forward = "%s/%s_R1_001.cor.fastq" %(newPath, name)
        reverse = "%s/%s_R2_001.cor.fastq" %(newPath, name)
        # There's an option to continue from checkpoints if the assembly is terminated prematurely within SPAdes,
        #  but the output directory must exist - if this directory exists, --continue, else don't --continue
        # /home/blais/Bioinformatics/SPAdes-3.1.1-Linux/bin/
        if os.path.isdir("%s/spades_output" % name):
            spadesRun = "spades.py -k 21,33,55,77,99,127 " \
                        "--careful --continue --only-assembler --pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
        else:
            spadesRun = "spades.py -k 21,33,55,77,99,127 --careful " \
                        "--only-assembler --pe1-1 %s --pe1-2 %s -o %s/spades_output 1>/dev/null" % (forward, reverse, newPath)
        # Run the command - subprocess.call would not run this command properly - no idea why - so using os.system instead
        # added 1>/dev/null to keep terminal output from being printed to screen
        os.system(spadesRun)
        # Print dots as per usual
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def contigFileFormatter(correctedFiles, path, metadata):
    """Changes the name of each contig from ">NODE_XXX_length..." to the name of the file ">OLC795_XXX..." """
    for name in correctedFiles:
        newPath = path + "/" + name
        # Ensures that the contigs file is present, but the renamed, manipulated file is not
        # os.remove("%s/%s_filteredAssembled.fasta" % (newPath, name))
        if os.path.isfile("%s/spades_output/contigs.fasta" % newPath) and not os.path.isfile("%s/%s_filteredAssembled.fasta" % (newPath, name)):
            # http://biopython.org/wiki/SeqIO#Input.2FOutput_Example_-_Filtering_by_sequence_length
            over200bp = []
            lengthCov = 0
            for record in SeqIO.parse(open("%s/spades_output/contigs.fasta" % newPath, "rU"), "fasta"):
                # Include only contigs greater than 200 bp in length
                if len(record.seq) >= 200:
                    # Add this record to our list
                    newID = re.sub("NODE", name, record.id)
                    lengthCov += (float(newID.split("_")[-5]) * float(newID.split("_")[-3]))
                    record.id = newID
                    record.name = ''
                    record.description = ''
                    over200bp.append(record)
            metadata[name]["2.Assembly"]["totalBasesxCoverage"] = lengthCov
            # print metadata[name]["2.Assembly"]["totalBasesxCoverage"]
            # print "Found %s long sequences" % len(over200bp)
            fileName = "%s/%s_filteredAssembled.fasta" % (newPath, name)
            formatted = open(fileName, "wb")
            SeqIO.write(over200bp, formatted, "fasta")
            formatted.close()
            # Move the files to a BestAssemblies folder
            assemblyPath = "%s/BestAssemblies" % path
            make_path(assemblyPath)
            fileName = "%s/%s_filteredAssembled.fasta" % (newPath, name)
            if not os.path.isfile("%s/%s" % (assemblyPath, fileName)):
                shutil.copy(fileName, assemblyPath)
    return metadata


def functionsGoNOW(correctedFiles, path, metadata):
    """Run the helper function"""
    print("\nAssembling reads.")
    spadesPrepProcesses(correctedFiles, path)
    updatedMetadata = contigFileFormatter(correctedFiles, path, metadata)
    return updatedMetadata