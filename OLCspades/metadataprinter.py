#!/usr/bin/env python
__author__ = 'adamkoziol'


class MetadataPrinter(object):

    def printmetadata(self):
        import json
        # Iterate through each sample in the analysis
        for sample in self.metadata.runmetadata.samples:
            # Open the metadata file to write
            with open('{}/{}_metadata.json'.format(sample.general.outputdirectory, sample.name), 'wb') as metadatafile:
                # Write the json dump of the object dump to the metadata file
                metadatafile.write(json.dumps(sample.dump(), sort_keys=True, indent=4, separators=(',', ': ')))

    def __init__(self, inputobject):
        self.metadata = inputobject
        self.printmetadata()
