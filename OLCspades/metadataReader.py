#!/usr/bin/env python
__author__ = 'adamkoziol'


class MetadataReader(object):

    def reader(self):
        import os
        import json
        from accessoryFunctions import GenObject, MetadataObject
        for sample in self.metadata:
            metadatafile = '{}{}/{}_metadata.json'.format(self.path, sample.name, sample.name)
            if os.path.isfile(metadatafile):
                with open(metadatafile) as metadatareport:
                    jsondata = json.load(metadatareport)
                # Create the metadata objects
                metadata = MetadataObject()
                # Set the name
                metadata.name = sample.name
                # Initialise the metadata categories as GenObjects created using the appropriate key
                metadata.run = GenObject(jsondata['run'])
                metadata.general = GenObject(jsondata['general'])
                metadata.commands = GenObject(jsondata['commands'])
                self.samples.append(metadata)

    def __init__(self, inputobject):
        # self.metadata = inputobject.runmetadata.samples
        self.metadata = inputobject.samples
        self.path = inputobject.path
        self.samples = []
        self.reader()
