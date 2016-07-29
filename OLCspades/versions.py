#!/usr/bin/env python
import sys

import Bio

from accessoryFunctions import *
from bowtie import *

__author__ = 'adamkoziol'


class Versions(object):
    def versions(self):
        for sample in self.metadata:
            # Initialise the attribute
            sample.software = GenObject()
            # Populate the versions of the software used
            ss = sample.software
            ss.python = self.python
            ss.arch = self.arch
            ss.blast = self.blast
            ss.bowtie2 = self.bowversion
            ss.samtools = self.samversion
            ss.qualimap = self.qualimap
            ss.mash = self.mash
            ss.prodigal = self.prodigal
            ss.pipeline = self.commit
            ss.quast = self.quast
            ss.spades = self.spades
            ss.bbmap = self.bbmap
            ss.fastqc = self.fastqc
            ss.blc2fastq = self.bcl2fastq
            ss.perl = self.perl
            ss.biopython = self.biopython
            ss.java = self.java
            # ss.docker = self.docker

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.commit = inputobject.commit
        # Determine the versions of the software used
        printtime('Populating metadata', self.start)
        self.python = sys.version.replace('\n', '')
        self.arch = ", ".join(os.uname())
        self.blast = get_version(['blastn', '-version']).split('\n')[0].split()[1]
        self.spades = get_version(['spades.py', '-v']).split('\n')[0].split()[1]
        self.bowversion = Bowtie2CommandLine(version=True)()[0].split('\n')[0].split()[-1]
        self.samversion = get_version(['samtools', '--version']).split('\n')[0].split()[1]
        # Qualimap seems to have an Java warning message that doesn't necessarily show up on every system
        # Only capture the line that starts with 'Qualimap'
        qualimaplist = get_version(['qualimap', '--help']).split('\n')
        for line in qualimaplist:
            if 'QualiMap' in line:
                self.qualimap = line.split()[1]
        self.mash = get_version(['mash']).split('\n')[1].split()[2]
        self.prodigal = get_version(['prodigal', '-v']).split('\n')[1].split()[1]
        self.quast = get_version(['quast.py']).split('\n')[1].split()[1]
        self.bbmap = get_version(['bbmap.sh']).split('\n')[1].split()[1]
        self.fastqc = get_version(['fastqc', '--version']).split('\n')[0].split()[1]
        self.bcl2fastq = get_version(['configureBclToFastq.pl']).split('\n')[-2].split()[2].split('/')[4].split('-')[1]
        self.perl = get_version(['perl', '-v']).split('\n')[1].split('This is ')[1]
        self.biopython = Bio.__version__
        self.java = get_version(['java', '-showversion']).split('\n')[0].split()[2].replace('"', '')
        # self.docker = get_version(['docker', 'version']).split('\n')[1].split()[1]
        self.versions()
