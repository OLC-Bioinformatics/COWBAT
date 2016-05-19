#!/usr/bin/env python
import os
from threading import Thread

from accessoryFunctions import *

__author__ = 'adamkoziol'


class Vtyper(object):

    def vtyper(self):
        """Setup and create  threads for blastn and xml path"""
        # Create the threads for the BLAST analysis
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                threads = Thread(target=self.epcr, args=())
                threads.setDaemon(True)
                threads.start()

        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'stx' in sample.general.datastore:
                    setattr(sample, self.analysistype, GenObject())
                    # Get the primers ready
                    sample[self.analysistype].primers = '{}{}/vtx_subtyping_primers.txt'\
                        .format(self.reffilepath, self.analysistype)
                    # Make the output path
                    sample[self.analysistype].reportdir = '{}/{}/'.format(sample.general.outputdirectory,
                                                                          self.analysistype)
                    make_path(sample[self.analysistype].reportdir)
                    outfile = sample[self.analysistype].reportdir + sample.name
                    # # Create the directory/sample.name (no extension)-containing variable
                    # linkfile = '{}{}'.format(sample[self.analysistype].reportdir, sample.name)
                    # sample[self.analysistype].linkfile = linkfile
                    # # Link the contigs file to the report dir
                    # try:
                    #     os.symlink(sample.general.filteredfile, '{}.fasta'.format(linkfile))
                    # # Except os errors
                    # except OSError as exception:
                    #     # If the os error is anything but directory exists, then raise
                    #     if exception.errno != errno.EEXIST:
                    #         raise
                    sample.commands.famap = 'famap -b {}.famap {}.fasta'.format(outfile, sample.general.filenoext)
                    sample.commands.fahash = 'fahash -b {}.hash {}.famap'.format(outfile, outfile)
                    # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
                    # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
                    # -m 10000 (Set variability for STS size for lookup),
                    # -n 1 (Set max allowed mismatches per primer for lookup)
                    # -g 0 (Set max allowed indels per primer for lookup),
                    # -G (Print alignments in comments), -o {output file}
                    sample.commands.epcr = 're-PCR -S {}.hash -r + -m 10000 -n 1 -g 0 -G -q -o {}.txt {}'\
                        .format(outfile, outfile, sample[self.analysistype].primers)
                    sample[self.analysistype].resultsfile = '{}.txt'.format(outfile)
                    self.epcrqueue.put((sample, outfile))
        self.epcrqueue.join()
        self.epcrparse()

    def epcr(self):
        from subprocess import call
        while True:
            sample, linkfile = self.epcrqueue.get()
            if not os.path.isfile('{}.famap'.format(linkfile)):
                call(sample.commands.famap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            if not os.path.isfile('{}.hash'.format(linkfile)):
                call(sample.commands.fahash, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            if not os.path.isfile('{}.txt'.format(linkfile)):
                call(sample.commands.epcr, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            self.epcrqueue.task_done()

    def epcrparse(self):
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                if 'stx' in sample.general.datastore:
                    # Initialise count - this allows for the population of vtyperresults with unique values
                    uniquecount = 0
                    # This populates vtyperresults with the verotoxin subtypes
                    toxinlist = []
                    if os.path.isfile(sample[self.analysistype].resultsfile):
                        epcrresults = open(sample[self.analysistype].resultsfile, 'r')
                        for result in epcrresults:
                            # Only the lines without a # contain results
                            if "#" not in result:
                                uniquecount += 1
                                # Split on \t
                                data = result.split('\t')
                                # The subtyping primer pair is the first entry on lines with results
                                vttype = data[0].split('_')[0]
                                # Push the name of the primer pair - stripped of anything after a _ to the dictionary
                                if vttype not in toxinlist:
                                    toxinlist.append(vttype)
                    # Create a string of the entries in list1 joined with ";"
                    toxinstring = ";".join(toxinlist)
                    # Save the string to the metadata
                    sample[self.analysistype].toxinprofile = toxinstring
                else:
                    setattr(sample, self.analysistype, GenObject())
                    sample[self.analysistype].toxinprofile = 'NA'
            else:
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].toxinprofile = 'NA'

    def __init__(self, inputobject, analysistype):
        from Queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.analysistype = analysistype
        self.reffilepath = inputobject.reffilepath
        self.epcrqueue = Queue()
        self.vtyper()
