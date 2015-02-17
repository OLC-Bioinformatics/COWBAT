__author__ = 'blais'

import os, glob, time, shutil, errno, re

# Initialise variables
formattedDate = ""
folderofInterest = ""

miSeqPath = "/media/miseq/MiSeqOutput"
backupPath = "/media/nas/backup/MiSeq/MiSeqOutput"
analysisPath = "/media/nas/akoziol/WGS_Spades"

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def run():
    # Go to the appropriate folder
    os.chdir(miSeqPath)

    # Glob all the folders in the directory
    miSeqfolders = glob.glob("*")

    # Sort the folders - as folders are named based on the date, the most recent folder will be the last entry in the list
    sortedFolders = sorted(miSeqfolders)

    # Determine which is the folder of interest, and parse the date from its title
    # - this will be used in creating folders later on
    if os.path.isdir(sortedFolders[-1]):
        folderofInterest = sortedFolders[-1]
        print folderofInterest
        date = folderofInterest.split("_")[0]
        formattedDate = "20" + date[:2] + "-" + date[2:4] + "-" + date[4:6]

    sampleCount = 0
    with open("%s/%s/SampleSheet.csv" % (miSeqPath, folderofInterest)) as sampleSheet:
        for entry in sampleSheet:
            if re.search("Sample_ID", entry):
                for subline in sampleSheet:
                    sampleCount += 1

    print "There are %s samples in this run" % sampleCount

    # As the point of this script is to ensure that the sequencing run is completed, and then copy the directory
    # to our NAS, we must ensure that the fastq files have been created/
    GZfiles = glob.glob("%s/%s/Data/Intensities/BaseCalls/*.gz" % (miSeqPath, folderofInterest))

    print len(GZfiles)
    while len(GZfiles) < 2 * (sampleCount + 1):
        print "Waiting for the run to finish"
        print len(GZfiles)
        time.sleep(300)
        GZfiles = glob.glob("%s/%s/Data/Intensities/BaseCalls/*.gz" % (miSeqPath, folderofInterest))

    # if not os.path.isdir("%s/%s" % (backupPath, folderofInterest)):
    #     shutil.copytree("%s/%s" % (miSeqPath, folderofInterest), "%s/%s" % (backupPath, folderofInterest))

    for files in GZfiles:
        fileName = os.path.split(files)[-1]
        if not os.path.isfile("%s/%s/%s" % (backupPath, folderofInterest, fileName)):
            make_path("%s/%s" % (backupPath, folderofInterest))
            print "Copying %s" % fileName
            shutil.copy(files, "%s/%s/%s" % (backupPath, folderofInterest, fileName))

    # Copy over the necessary files to perform metadata analyses
    if not os.path.isfile("%s/%s/SampleSheet.csv" % (backupPath, folderofInterest)):
        shutil.copy("%s/%s/SampleSheet.csv" % (miSeqPath, folderofInterest),
                    "%s/%s/SampleSheet.csv" % (backupPath, folderofInterest))
        shutil.copy("%s/%s/RunInfo.xml" % (miSeqPath, folderofInterest),
                    "%s/%s/RunInfo.xml" % (backupPath, folderofInterest))
        if os.path.isfile("%s/%s/GenerateFASTQRunStatistics.xml" % (miSeqPath, folderofInterest)):
            shutil.copy("%s/%s/GenerateFASTQRunStatistics.xml" % (miSeqPath, folderofInterest),
                        "%s/%s/GenerateFASTQRunStatistics.xml" % (backupPath, folderofInterest))

    analysisFiles = glob.glob("%s/%s/*" % (backupPath, folderofInterest))

    # Make the path in the analysis folder
    make_path("%s/%s" % (analysisPath, formattedDate))

    for aFile in analysisFiles:
        if not re.search("Undetermined_S0_L001_", aFile):
            aFileName = os.path.split(aFile)[-1]
            folderName = re.split("_S\d+_L001", aFileName)[0]
            # print "%s/%s/%s" % (analysisPath, formattedDate, folderName)
            if not os.path.isdir("%s/%s/%s" % (analysisPath, formattedDate, folderName)) and re.search(".gz", aFileName)\
                    and not os.path.isfile("%s/%s/%s" % (analysisPath, formattedDate, aFileName)):
                os.symlink(aFile, "%s/%s/%s" % (analysisPath, formattedDate, aFileName))
            elif not os.path.isfile("%s/%s/%s" % (analysisPath, formattedDate, aFileName)) and not re.search(".gz", aFileName):
                os.symlink(aFile, "%s/%s/%s" % (analysisPath, formattedDate, aFileName))

    os.chdir("%s/%s" % (analysisPath, formattedDate))

# run()