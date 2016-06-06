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
                size = os.stat(metadatafile).st_size
                if size != 0:
                    try:
                        with open(metadatafile) as metadatareport:
                            jsondata = json.load(metadatareport)
                        # Create the metadata objects
                        metadata = MetadataObject()
                        # Initialise the metadata categories as GenObjects created using the appropriate key
                        for attr in jsondata:
                            if not isinstance(jsondata[attr], dict):
                                setattr(metadata, attr, jsondata[attr])
                            else:
                                setattr(metadata, attr, GenObject(jsondata[attr]))
                        # Set the name
                        metadata.name = sample.name
                        self.samples.append(metadata)
                    except ValueError:
                        self.samples.append(sample)
            else:
                self.samples.append(sample)

    def __init__(self, inputobject):
        self.metadata = inputobject.samples
        self.path = inputobject.path
        self.samples = []
        self.reader()
