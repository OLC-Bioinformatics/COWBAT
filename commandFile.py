__author__ = 'akoziol'

import os, json
from collections import defaultdict
import metadataFiller

def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


def commands(sampleNames, path, numReads, metadata):
    """Pulls all the commands performed in a previous run of the pipeline and enters them into a dictionary"""
    # Make an empty defaultdict dictionary
    performedCommands = defaultdict(make_dict)
    programVersions = defaultdict(make_dict)
    strReads = "paired-end"
    # Convert numReads (1 or 2) to a string (single or paired-end)
    if int(numReads) == 1:
        strReads = "single"
    # For every strain sequenced, search to see if the metadata json report has been created.
    for name in sampleNames:
        if os.path.isfile("%s/%s/%s_metadataCollection.json" % (path, name, name)):
            # Open the json file
            try:
                with open("%s/%s/%s_metadataReport.json" % (path, name, name)) as jsonReport:
                    # Load the data
                    jsonData = json.load(jsonReport)
                    # Iterate through the commands stored in the 7.Pipeline section
                    if "7.PipelineCommands" in jsonData:
                        for command in jsonData["7.PipelineCommands"]:
                            # If there is an "N/A" entry, then the command was either not performed, or not stored
                            # the way the pipeline will now be run is that if the command is not stored, then the
                            # section of the pipeline missing the command will be re-performed. Hopefully, this will
                            # allow for quickly (re-)performing necessary analysis.
                            if not "N/A" in jsonData["7.PipelineCommands"][command] and not "defaultdict" in jsonData["7.PipelineCommands"][command]:
                                # Populate performedCommands with the commands previously performed
                                performedCommands[name][str(command)] = str(jsonData["7.PipelineCommands"][command])
                    # Perform the same actions for the versions section of the json file
                    if "8.PipelineVersions" in jsonData:
                        for command in jsonData["8.PipelineVersions"]:
                            # If there is an "N/A" entry, then the command was either not performed, or not stored
                            # the way the pipeline will now be run is that if the command is not stored, then the
                            # section of the pipeline missing the command will be re-performed. Hopefully, this will
                            # allow for quickly (re-)performing necessary analysis.
                            if not "N/A" in jsonData["8.PipelineVersions"][command] and not "defaultdict" in jsonData["8.PipelineVersions"][command]:
                                # Populate performedCommands with the commands previously performed
                                programVersions[name][str(command)] = str(jsonData["8.PipelineVersions"][command])
                    if "3.Run" in jsonData:
                        if "PairedEnd/SingleReads" in jsonData["3.Run"]:
                            if not "N/A" in jsonData["3.Run"]["PairedEnd/SingleReads"] and not "defaultdict" in jsonData["3.Run"]["PairedEnd/SingleReads"]:
                                strReads = jsonData["3.Run"]["PairedEnd/SingleReads"]
                                metadata[name]["3.Run"]["PairedEnd/SingleReads"] = strReads
                        else:
                            metadata[name]["3.Run"]["PairedEnd/SingleReads"] = strReads
            except (RuntimeError, TypeError, NameError, IOError):
                with open("%s/%s/%s_metadataCollection.json" % (path, name, name)) as jsonReport:
                    # Load the data
                    jsonData = json.load(jsonReport)
                    # Iterate through the commands stored in the 7.Pipeline section
                    if "7.PipelineCommands" in jsonData:
                        for command in jsonData["7.PipelineCommands"]:
                            # If there is an "N/A" entry, then the command was either not performed, or not stored
                            # the way the pipeline will now be run is that if the command is not stored, then the
                            # section of the pipeline missing the command will be re-performed. Hopefully, this will
                            # allow for quickly (re-)performing necessary analysis.
                            if not "N/A" in jsonData["7.PipelineCommands"][command] and not "defaultdict" in jsonData["7.PipelineCommands"][command]:
                                # Populate performedCommands with the commands previously performed
                                performedCommands[name][str(command)] = str(jsonData["7.PipelineCommands"][command])
                    if "8.PipelineVersions" in jsonData:
                        for command in jsonData["8.PipelineVersions"]:
                            # If there is an "N/A" entry, then the command was either not performed, or not stored
                            # the way the pipeline will now be run is that if the command is not stored, then the
                            # section of the pipeline missing the command will be re-performed. Hopefully, this will
                            # allow for quickly (re-)performing necessary analysis.
                            if not "N/A" in jsonData["8.PipelineVersions"][command] and not "defaultdict" in jsonData["8.PipelineVersions"][command]:
                                # Populate performedCommands with the commands previously performed
                                programVersions[name][str(command)] = str(jsonData["8.PipelineVersions"][command])
                    if "3.Run" in jsonData:
                        if "PairedEnd/SingleReads" in jsonData["3.Run"]:
                            if not "N/A" in jsonData["3.Run"]["PairedEnd/SingleReads"] and not "defaultdict" in jsonData["3.Run"]["PairedEnd/SingleReads"]:
                                strReads = jsonData["3.Run"]["PairedEnd/SingleReads"]
                                metadata[name]["3.Run"]["PairedEnd/SingleReads"] = strReads
                        else:
                            metadata[name]["3.Run"]["PairedEnd/SingleReads"] = strReads
        else:
            metadata[name]["3.Run"]["PairedEnd/SingleReads"] = strReads
    # Return the variables
    return performedCommands, programVersions, strReads, metadata
