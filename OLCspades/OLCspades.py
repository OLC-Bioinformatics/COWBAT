#!/usr/bin/env python
import os
import runMetadata
import offhours
import json
__author__ = 'adamkoziol'


class KeyboardInterruptError(Exception):
    pass


class RunSpades(object):
    def __init__(self, args):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.path = os.path.join(args['path'], "")
        self.numreads = args['numReads']
        self.offhours = args['offHours']
        self.reffilepath = os.path.join(args['referenceFilePath'], "")
        self.kmers = args['kmerRange']
        self.customsamplesheet = args['customSampleSheet']
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        if args['threads']:
            self.cpus = args['threads']
        else:
            self.cpus = os.popen("awk '/^processor/ { N++} END { print N }' /proc/cpuinfo").read().rstrip()
        # Assertions to ensure that the provided variables are valid
        assert os.path.isdir(self.path), u'Output location is not a valid directory {0!r:s}'.format(self.path)
        assert os.path.isdir(self.reffilepath), u'Output location is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        # Initialise the metadata object
        self.runmetadata = ""

    def assembly(self):
        """Performs the assembly"""
        print("Extracting metadata from sequencing run.")
        # Populate the runmetadata object by parsing the SampleSheet.csv, GenerateFASTQRunStatistics.xml, and
        # RunInfo.xml files
        self.runmetadata = runMetadata.Metadata(self)
        # Run the offhours fastq creation script
        if self.offhours:
            offhours.Offhours(self)
        # Extract the flowcell ID and the instrument name if the RunInfo.xml file was provided
        self.runmetadata.parseruninfo()
        # print json.dumps([x.dump() for x in self.runmetadata.samples], sort_keys=True, indent=4, separators=(',', ': '))


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Get the current commit of the pipeline from git
    blerg = 'blarg'
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Assemble genomes from Illumina fastq files')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s commit {}'.format(blerg))
    parser.add_argument('path',  help='Specify path')
    parser.add_argument('-n', metavar='numreads', default=2, help='Specify the number of reads. Paired-reads: 2, '
                        'unpaired-reads: 1. Default is paired-end')
    parser.add_argument('-t', metavar='threads', help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-o', '--offHours', action='store_true', help='Optionally run the off-hours module that will '
                        'search for MiSeq runs in progress, wait until the run is complete, and assemble the run')
    parser.add_argument('-F', '--FastqCreation', action='store_true', help='Optionally run the fastq creation module'
                        'that will search for MiSeq runs in progress, run bcl2fastq to create fastq files, and '
                        'assemble the run')
    parser.add_argument('-m', metavar='miSeqPath', help='Path of the folder containing MiSeq run data folder')
    parser.add_argument('-f', metavar='miseqfolder', help='Name of the folder containing MiSeq run data')
    parser.add_argument('-r1', metavar='readLengthForward', help='Length of forward reads to use. Can specify'
                        '"full" to take the full length of forward reads specified on the SampleSheet')
    parser.add_argument('-r2', metavar='readLengthReverse', default=0, help='Length of reverse reads to use. '
                        'Can specify "full" to take the full length of reverse reads specified on the SampleSheet')
    parser.add_argument('-r', metavar='referenceFilePath', default="/spades_pipeline/SPAdesPipelineFiles",
                        help='Provide the location of the folder containing the pipeline accessory files '
                        '(reference genomes, MLST data, etc.')
    parser.add_argument('-k', metavar='kmerRange', default='21,33,55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 21,33,55,77,99,127')
    parser.add_argument('-c', metavar='customSampleSheet', help='Path of folder containing a custom sample '
                        'sheet (still must be named "SampleSheet.csv")')

    # Get the arguments into a list
    arguments = vars(parser.parse_args())

    # Run the pipeline
    output = RunSpades(arguments)
    output.assembly()
