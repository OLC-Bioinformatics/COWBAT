#!/usr/bin/env python
import json
__author__ = 'adamkoziol'


class MetadataPrinter(object):

    def printmetadata(self):
        # Iterate through each sample in the analysis
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Set the name of the json file
                jsonfile = '{}/{}_metadata.json'.format(sample.general.outputdirectory, sample.name)
                try:
                    # Open the metadata file to write
                    with open(jsonfile, 'wb') as metadatafile:
                        # Write the json dump of the object dump to the metadata file
                        json.dump(sample.dump(), metadatafile, sort_keys=True, indent=4, separators=(',', ': '))
                except IOError:
                    # print json.dumps(sample.datastore, sort_keys=True, indent=4, separators=(',', ': '))
                    print sample.datastore
                    raise

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.printmetadata()
