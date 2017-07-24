#!/usr/bin/env python
from accessoryFunctions import *
from glob import glob
__author__ = 'adamkoziol'


class ObjectCreation(object):

    def createobject(self):
        printtime('Finding sequence files', self.start)
        # Find all the .fastq files in the sequence path
        filelist = glob('{}*.fastq*'.format(self.sequencepath))
        if filelist:
            self.extension = '.fastq'
            self.filehandler(filelist)
        else:
            filelist = glob('{}*.fa*'.format(self.sequencepath))
            self.extension = '.fasta'
            self.filehandler(filelist)

    def filehandler(self, filelist):
        # Extract the base name of the globbed name + path provided
        if self.extension == '.fastq':
            names = map(lambda x: os.path.split(x)[1], filer(filelist))
        else:
            names = map(lambda x: os.path.split(x)[1].split('.')[0], filelist)
        # Iterate through the names of the fastq files
        for name in sorted(names):
            # Set the name
            metadata = MetadataObject()
            metadata.name = name
            # Set the destination folder
            outputdir = os.path.join(self.sequencepath, name)
            # Make the destination folder
            make_path(outputdir)
            # Get the fastq files specific to the fastqname
            specific = glob('{}{}*{}*'.format(self.sequencepath, name, self.extension))
            # Link the files to the output folder
            try:
                # Link the .gz files to :self.path/:filename
                map(lambda x: os.symlink(x, '{}/{}'.format(outputdir, os.path.split(x)[1])), specific)
            # Except os errors
            except OSError as exception:
                # If there is an exception other than the file exists, raise it
                if exception.errno != errno.EEXIST:
                    raise
            # Initialise the general and run categories
            metadata.general = GenObject()
            metadata.commands = GenObject()
            # Populate the .fastqfiles category of :self.metadata
            metadata.general.fastqfiles = [fastq for fastq in glob('{}/{}*{}*'.format(outputdir, name, self.extension))
                                           if 'trimmed' not in fastq]
            # Add the output directory to the metadata
            metadata.general.outputdirectory = outputdir

            # Find the data files corresponding to the sample
            datafiles = glob('{}{}*.csv'.format(self.datapath, metadata.name))
            # Assign attributes to the files depending on whether they are abundance files or not
            for datafile in datafiles:
                if 'abundance' in datafile:
                    metadata.general.abundancefile = datafile
                else:
                    metadata.general.assignmentfile = datafile

            # Append the metadata to the list of samples
            self.samples.append(metadata)

    def __init__(self, inputobject):
        self.samples = list()
        self.path = inputobject.path
        self.sequencepath = inputobject.sequencepath
        try:
            self.datapath = inputobject.datapath
        except AttributeError:
            self.datapath = False
        if self.datapath:
            assert os.path.isdir(self.datapath), u'Data location supplied is not a valid directory {0!r:s}' \
                .format(self.datapath)
        self.start = inputobject.start
        self.extension = ''
        self.createobject()
