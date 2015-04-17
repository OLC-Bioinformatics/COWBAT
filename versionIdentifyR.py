__author__ = 'akoziol'

import os, subprocess


def pipelineMetadata(path, metadata, sampleNames):
    # Update these values as required
    #Quast
    quastOutput = subprocess.Popen(["quast.py"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = quastOutput.communicate()
    quastVer = err.split("\n")[1].strip("Version ").replace(" ,", ",").replace(",", ";")

    # SMALT
    smaltOutput = subprocess.Popen(["smalt"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = smaltOutput.communicate()
    smaltVer = out.split("\n")[2].strip(" ").strip("(").strip(")").strip("version: ")

    # Samtools
    samOutput = subprocess.Popen(["samtools"], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, shell=True)
    out, err = samOutput.communicate()
    samVer = err.split("\n")[2].strip("Version: ")

    # Blast
    blastnOutput = subprocess.Popen(["blastn -version"], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
    out, err = blastnOutput.communicate()
    blastVer = out.split("\n")[0].strip("blastn: ")

    # Pipeline commits
    spadesLocation = subprocess.Popen(["which SPAdesPipeline.py"], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell=True)
    out, err = spadesLocation.communicate()
    spadesLoc = out
    spadesFolder = os.path.dirname(spadesLoc)
    spadesCommitOut = subprocess.Popen(["git --git-dir %s/.git log --date short --pretty=format:'%%h - %%an - %%ad' -1"
                                        % spadesFolder], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = spadesCommitOut.communicate()
    spadesCommit = out

    for name in sampleNames:

        with open("%s/%s/spades_output/spades.log" % (path, name)) as spadesLog:
            for line in spadesLog:
                if "SPAdes version" in line:
                    spadesVer = line.strip("SPAdes version: ").rstrip()
                elif "Python version" in line:
                    pythonVer = line.strip("Python version: ").rstrip()
                elif "OS" in line:
                    OS = line.strip("OS: ").rstrip()

        quakeVer = "0.3"
        # commit = "b737e2c52f59c541062a5f71c5e08087ee238c1e 2015-02-17 Adam Koziol"

        if not os.path.isfile("%s/%s/%s_programVersions.tsv" % (path, name, name)):

            versions = open("%s/%s/%s_programVersions.tsv" % (path, name, name), "wb")
            versions.write("spadesVersion\tquastVersion\tquakeVersion\tSmaltVersion\tSamtools\tBlastVersion\t"
                       "PythonVersion\tOS\tCommit\n")
            versions.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (spadesVer, quastVer, quakeVer, smaltVer, samVer, blastVer, pythonVer, OS, spadesCommit))

            # print name, quastVer, smaltVer, samVer, blastVer, spadesVer, pythonVer, OS
            metadata[name]["8.PipelineVersions"]["SPAdesVersion"] = spadesVer
            metadata[name]["8.PipelineVersions"]["QUASTVersion"] = quastVer
            metadata[name]["8.PipelineVersions"]["QuakeVersion"] = quakeVer
            metadata[name]["8.PipelineVersions"]["SmaltVersion"] = smaltVer
            metadata[name]["8.PipelineVersions"]["SamtoolsVersion"] = samVer
            metadata[name]["8.PipelineVersions"]["BlastVersion"] = blastVer
            metadata[name]["8.PipelineVersions"]["PythonVersion"] = pythonVer
            metadata[name]["8.PipelineVersions"]["OS"] = OS
            metadata[name]["8.PipelineVersions"]["PipelineVersion"] = spadesCommit
        else:
            pipelineVersions = open("%s/%s/%s_programVersions.tsv" % (path, name, name)).readlines()[1]
            versions = pipelineVersions.split("\t")
            metadata[name]["8.PipelineVersions"]["SPAdesVersion"] = versions[0]
            metadata[name]["8.PipelineVersions"]["QUASTVersion"] = versions[1]
            metadata[name]["8.PipelineVersions"]["QuakeVersion"] = versions[2]
            metadata[name]["8.PipelineVersions"]["SmaltVersion"] = versions[3]
            metadata[name]["8.PipelineVersions"]["SamtoolsVersion"] = versions[4]
            metadata[name]["8.PipelineVersions"]["BlastVersion"] = versions[5]
            metadata[name]["8.PipelineVersions"]["PythonVersion"] = versions[6]
            metadata[name]["8.PipelineVersions"]["OS"] = versions[7]
            metadata[name]["8.PipelineVersions"]["PipelineVersion"] = versions[8].rstrip()

        if not metadata[name]["4.Correction"]["ForwardValidatedReads"]:
            metadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardValidatedReads"]:
            metadata[name]["4.Correction"]["ForwardValidatedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardCorrectedReads"]:
            metadata[name]["4.Correction"]["ForwardCorrectedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardTrimmedReads"]:
            metadata[name]["4.Correction"]["ForwardTrimmedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"]:
            metadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"]:
            metadata[name]["4.Correction"]["ForwardRemovedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ForwardTrimmedOnlyReads"]:
            metadata[name]["4.Correction"]["ReverseValidatedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ReverseCorrectedReads"]:
            metadata[name]["4.Correction"]["ReverseCorrectedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ReverseTrimmedReads"]:
            metadata[name]["4.Correction"]["ReverseTrimmedReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ReverseTrimmedOnlyReads"]:
            metadata[name]["4.Correction"]["ReverseTrimmedOnlyReads"] = "N/A"
        if not metadata[name]["4.Correction"]["ReverseRemovedReads"]:
            metadata[name]["4.Correction"]["ReverseRemovedReads"] = "N/A"
    return metadata
