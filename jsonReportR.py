__author__ = 'akoziol'

import json, shutil, os, errno


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def jsonR(sampleNames, path, metadata, fileType):
    reportPath = "%s/reports" % path
    make_path(reportPath)
    for name in sampleNames:
            newPath = path + "/" + name
            reportName = "%s_metadata%s.json" % (name, fileType)
            JSONreport = open("%s/%s" % (newPath, reportName), "wb")
            # Print the JSON data to file
            output = json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
            JSONreport.write(output)
            JSONreport.close()
            # Move all the reports to a common directory
            shutil.copyfile("%s/%s" % (newPath, reportName), "%s/%s" % (reportPath, reportName))