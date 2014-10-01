__author__ = 'akoziol'

import json
import errno
import os
import shutil
import re


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def functionsGoNOW(sampleNames, metadata, path):
    """As outputting a JSON file is so straightforward, helper functions were unnecessary"""
    print "\nCreating reports."
    for name in sampleNames:
        newPath = path + "/" + name
        reportName = "%s_metadataReport.json" % name
        JSONreport = open("%s/%s" % (newPath, reportName), "wb")
        # Print the JSON data to file
        output = json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
        JSONreport.write(output)
        JSONreport.close()
        # Move all the reports to a common directory
        reportPath = "%s/reports" % path
        make_path(reportPath)
        shutil.copy("%s/%s" % (newPath, reportName), reportPath)
