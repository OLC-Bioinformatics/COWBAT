__author__ = 'akoziol'

import json
import re

def functionsGoNOW(sampleNames, metadata, path):
    print "Creating reports."
    for name in sampleNames:
        newPath = path + "/" + name
        JSONreport = open("%s/%s_metadataReport.json" % (newPath, name), "wb")
        output = json.dumps(metadata[name], sort_keys=True, indent=4)
        JSONreport.write(output)
        JSONreport.close()
    # print path
