#!/usr/bin/env python
from glob import glob
from subprocess import call
from threading import Lock

from Bio.Sequencing.Applications import SamtoolsViewCommandline, SamtoolsSortCommandline

from accessoryFunctions import *
from bowtie import *

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

threadlock = Lock()
__author__ = 'mike knowles, adamkoziol'


class QualiMap(object):

    def bowtie(self):
        from threading import Thread
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.align, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # main = lambda (x, y): (y, ",".join(getattr(sample.general, x))) if hasattr(sample.general, x) else None
            sagen = sample.general
            if sagen.bestassemblyfile != "NA":
                sagen.QualimapResults = '{}/qualimap_results'.format(sagen.outputdirectory)
                # Set the results folder
                # Create this results folder if necessary
                make_path(sagen.QualimapResults)
                # sagen.bamfile = sample.general.filenoext + '_sorted'
                sagen.sortedbam = '{}/{}_sorted.bam'.format(sagen.QualimapResults, sample.name)
                filenoext = sagen.filteredfile.split('.')[0]
                sagen.filenoext = filenoext
                sagen.bowtie2results = os.path.join(sagen.QualimapResults, sample.name)
                # Use fancy new bowtie2 wrapper
                bowtie2build = Bowtie2BuildCommandLine(reference=sagen.bestassemblyfile,
                                                       bt2=sagen.bowtie2results)
                sample.mapping.BamFile = sagen.bowtie2results + "_sorted.bam"
                # SAMtools sort v1.3 has different run parameters
                if self.samversion < "1.3":
                    samsort = SamtoolsSortCommandline(input_bam="-", out_prefix=sample.mapping.BamFile)
                else:
                    samsort = SamtoolsSortCommandline(input_bam=sample.mapping.BamFile,
                                                      o=True,
                                                      out_prefix="-")
                samtools = [SamtoolsViewCommandline(b=True, S=True, input_file="-"), samsort]
                indict = {'D': 5, 'R': 1, 'num_mismatches': 0, 'seed_length': 22, 'i_func': "S,0,2.50"}
                #  and sample.run.forwardlength > 50
                if len(sample.general.assemblyfastq) == 2:
                    indict.update({'m1': sample.general.assemblyfastq[0], 'm2': sample.general.assemblyfastq[1]})
                else:
                    indict.update({'U': sample.general.assemblyfastq[0]})
                bowtie2align = Bowtie2CommandLine(bt2=sagen.bowtie2results,
                                                  threads=self.cpus,
                                                  samtools=samtools,
                                                  **indict)

                sample.commands.bowtie2align = str(bowtie2align)
                sample.commands.bowtie2build = str(bowtie2build)
                self.bowqueue.put((sample, bowtie2build, bowtie2align))
            else:
                sample.commands.samtools = "NA"
        self.bowqueue.join()

    def align(self):
        while True:
            sample, bowtie2build, bowtie2align = self.bowqueue.get()
            if sample.general.bestassemblyfile != 'NA':
                if not os.path.isfile(sample.mapping.BamFile) and not os.path.isfile(sample.mapping.BamFile + ".bz2"):
                    stdout = StringIO()
                    for func in bowtie2build, bowtie2align:
                        stdout.close()
                        # Use cStringIO streams to handle bowtie output
                        stdout, stderr = map(StringIO, func(cwd=sample.general.QualimapResults))
                        if stderr:
                            # Write the standard error to log, bowtie2 puts alignment summary here
                            with open(os.path.join(sample.general.QualimapResults, "bowtieSamtools.log"), "ab+") as log:
                                log.writelines(logstr(func, stderr.getvalue(), stdout.getvalue()))
                        stderr.close()
                        # stdout will be the SAM file from alignment
            # For different alignment
            sam = sample.general.bowtie2results + ".sam"
            if os.path.isfile(sam):
                # PIPE stdout to stdin of samtools view then sort (only outputing sorted bam)
                # SAMtools sort v1.3 has different run parameters
                if self.samversion < "1.3":
                    samsort = SamtoolsSortCommandline(input_bam="-", out_prefix=sample.mapping.BamFile[:-4])
                else:
                    samsort = SamtoolsSortCommandline(input_bam=sample.mapping.BamFile, o=True, out_prefix="-")
                # Use cStringIO streams to handle bowtie output
                stdout = StringIO()
                for func in [SamtoolsViewCommandline(b=True, S=True, input_file=sample.mapping.BamFile), samsort]:
                    # Use closing contextmanager for handle __exit__() as close()
                    stdout, stderr = map(StringIO, func(stdin=stdout.getvalue()))
                    # Write the standard error to log
                    with open(os.path.join(sample.general.QualimapResults, "samtools.log"), "ab+") as log:
                        log.writelines(logstr(func, stderr.getvalue()))
                    stderr.close()
                stdout.close()
            self.mapper(sample)
            # Signal to the queue that the job is done
            self.bowqueue.task_done()

    def mapper(self, sample):
        if sample.general.bestassemblyfile != "NA":
            # Define the Qualimap log and report files
            reportfile = os.path.join(sample.general.QualimapResults, 'genome_results.txt')
            # Define the Qualimap call
            qualimapcall = 'qualimap bamqc -bam {} -outdir {}'.format(sample.general.sortedbam,
                                                                      sample.general.QualimapResults)
            sample.commands.qualimap = qualimapcall
            # Initialise a dictionary to hold the Qualimap results
            qdict = dict()
            # If the report file doesn't exist, run Qualimap, and print logs to the log file
            if not os.path.isfile(reportfile):
                call(sample.commands.qualimap, shell=True, stdout=self.fnull, stderr=self.fnull)
            try:
                with open(reportfile) as report:
                    # Read the report
                    for line in report:
                        # Sanitise the keys and values using self.analyze
                        key, value = self.analyze(line)
                        # If the keys and values exist, enter them into the dictionary
                        if (key, value) != (None, None):
                            qdict[key] = value
            except IOError:
                self.bowqueue.task_done()
                self.bowqueue.join()
                raise

            # If there are values in the dictionary
            if qdict:
                # Make new category for Qualimap results and populate this category with the report data
                setattr(sample, "mapping", GenObject(qdict))

    @staticmethod
    def analyze(line):
        # Split on ' = '
        if ' = ' in line:
            key, value = line.split(' = ')
            # Replace occurrences of
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Should we keep comma separation?
            value = value.replace(",", "").replace(" ", "").rstrip()
        # Otherwise set the keys and values to None
        else:
            key, value = None, None
        return key, value

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        # Define /dev/null
        self.fnull = open(os.devnull, 'wb')
        self.samversion = get_version(['samtools']).split('\n')[2].split()[1]
        # Initialise queues
        self.mapqueue = Queue(maxsize=self.cpus)
        self.qqueue = Queue(maxsize=self.cpus)
        self.bowqueue = Queue(maxsize=self.cpus)
        printtime('Aligning reads with bowtie2 for Qualimap', self.start)
        self.bowtie()


if __name__ == '__main__':
    class Parser(object):

        def associate(self):
            from accessoryFunctions import GenObject, MetadataObject
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = [fasta for fasta in sorted(glob('{}*.fa*'.format(self.assemblypath)))
                            if '.fastq' not in fasta]
            for strain in self.strains:
                # Extract the name of the strain from the path and file extension
                strainname = os.path.split(strain)[1].split('.')[0]
                # Find the corresponding fastq files for each strain
                fastq = sorted(glob('{}{}*fastq*'.format(self.fastqpath, strainname)))
                # Ensure that fastq files are present for each assembly
                assert fastq, 'Cannot find fastq files for strain {}'.format(strainname)
                # Create the object
                metadata = MetadataObject()
                # Set the .name attribute to be the file name
                metadata.name = strainname
                # Create the .general attribute
                metadata.general = GenObject()
                # Set the .general.filteredfile file to be the name and path of the sequence file
                # metadata.general.filteredfile = strain
                # Set the path of the assembly file
                metadata.general.bestassembliespath = self.assemblypath
                # Populate the .fastqfiles category of :self.metadata
                metadata.general.trimmedfastqfiles = fastq
                # Create the output directory path
                metadata.general.outputdirectory = '{}{}'.format(self.path, strainname)
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            from argparse import ArgumentParser
            import multiprocessing
            parser = ArgumentParser(description='Calculates coverage depth by mapping FASTQ reads against assemblies')
            parser.add_argument('-p', '--path',
                                default=os.getcwd(),
                                help='Specify the path of the folder that either contains the files of interest, or'
                                     'will be used to store the outputs')
            parser.add_argument('-a', '--assemblies',
                                help='Path to a folder of assemblies. If not provided, the script will look for .fa'
                                     'or .fasta files in the path')
            parser.add_argument('-f', '--fastq',
                                help='Path to a folder of fastq files. If not provided, the script will look for '
                                     'fastq or .fastq.gz files in the path')
            parser.add_argument('-t', '--threads',
                                help='Number of threads. Default is the number of cores in the system')
            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.assemblypath = os.path.join(args.assemblies, '') if args.assemblies else self.path
            self.fastqpath = os.path.join(args.fastq, '') if args.fastq else self.path
            # Use the argument for the number of threads to use, or default to the number of cpus in the system
            self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
            # Initialise variables
            self.strains = []
            self.samples = []

            # Associate the assemblies and fastq files in a metadata object
            self.associate()

    class MetadataInit(object):
        def __init__(self, start):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.assemblypath = self.runmetadata.assemblypath
            self.fastqpath = self.runmetadata.fastqpath
            self.starttime = start
            self.cpus = self.runmetadata.cpus
            # Run the analyses - the extra set of parentheses is due to using the __call__ method in the class
            QualiMap(self)

    # Run the class
    from time import time
    starttime = time()
    MetadataInit(starttime)
    printtime('Assembly and characterisation complete', starttime)
