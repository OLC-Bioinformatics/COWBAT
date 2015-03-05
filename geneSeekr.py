__author__ = 'akoziol'

from multiprocessing import Pool
from Bio.Blast.Applications import NcbiblastnCommandline
from collections import defaultdict
from Bio.Blast import NCBIXML
import errno, os, glob, json, sys, shutil

# Import ElementTree - try first to import the faster C version, if that doesn't
# work, try to import the regular version
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
geneSeekrResults = defaultdict(make_dict)
qualityResults = defaultdict(make_dict)
vtyperResults = defaultdict(make_dict)
profileData = defaultdict(make_dict)
MLSTresults = defaultdict(make_dict)
MLSTseqType = defaultdict(make_dict)

# Count variable is used in the dotter function
count = 0


def dotter():
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n.')
        count = 0


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def univecScreen(name, path, refFilePath):
    """Runs all the sequences through the Illumina-specific UniVec database - output is a .tsv file
    in the path/name/univecscreen folder"""
    # Set the primer, query, and out variables for the BLAST search
    primers = "%s/illuminaPrimers/IlluminaPrimers.fa" % refFilePath
    out = "%s/%s/univecscreen/%s.tsv" % (path, name, name)
    query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
    # Make the path (if necessary)
    make_path("%s/%s/univecscreen" % (path, name))
    # If the analysis hasn't already been performed
    if not os.path.isfile(out):
        # I'm running this with an evalue of 0.1 - the default VecScreen parameters call for an evalue of 700,
        # but since this is a reduced database, I reduced the evalue as well. I'm also using -task blastn-short
        # instead of blastn, as all the probes are very short. Additionally, as it was causing blast to return
        # no results, I removed -searchsp 1750000000000 from the call
        blastn = NcbiblastnCommandline(task="blastn-short",
                                       query=query,
                                       db=primers,
                                       reward=1,
                                       penalty=-5,
                                       gapopen=3,
                                       gapextend=3,
                                       dust="yes",
                                       soft_masking="true",
                                       evalue=0.1,
                                       outfmt=6,
                                       out=out)
        stdout, stderr = blastn()


def vtyper(name, path, refFilePath):
    """Uses ePCR to subtype genomic verotoxin sequences"""
    metadata = defaultdict(make_dict)
    # Because this is a multi-processed function, getting this dictionary to be properly formatted was difficult
    # in the end, I used if output, return output statements in the appropriate functions
    # global vtyperResults
    # Initialise count - this allows for the population of vtyperResults with unique values
    count = 0
    # Get the primers ready
    primers = "%s/vtyper/vtx_subtyping_primers.txt" % refFilePath
    # Make the output path
    vtyperPath = "%s/%s/vtyper" % (path, name)
    make_path(vtyperPath)
    # Copy the contigs file to the output path for ease of ePCR indexing
    shutil.copy("%s/%s/%s_filteredAssembled.fasta" % (path, name, name), "%s/%s/vtyper/" % (path, name))
    # For ease of processing set the contigs file as a variable
    query = "%s/%s_filteredAssembled.fasta" % (vtyperPath, name)
    # File checks to save re-performing processes
    if not os.path.isfile("%s/%s.famap" % (vtyperPath, name)):
        # famap and fahash are pre-processing steps performed on the query file
        famap = "famap -b %s/%s.famap %s 2>/dev/null" % (vtyperPath, name, query)
        os.system(famap)
    if not os.path.isfile("%s/%s.hash" % (vtyperPath, name)):
        fahash = "fahash -b %s/%s.hash %s/%s.famap 2>/dev/null" % (vtyperPath, name, vtyperPath, name)
        os.system(fahash)
    if not os.path.isfile("%s/%s.txt" % (vtyperPath, name)):
        # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
        # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
        # -m 10000 (Set variability for STS size for lookup), -n 1 (Set max allowed mismatches per primer for lookup)
        # -g 0 (Set max allowed indels per primer for lookup), -G (Print alignments in comments), -o {output file}
        ePCR = "re-PCR -S %s/%s.hash -r + -m 10000 -n 1 -g 0 -G -q -o %s/%s.txt %s 2>/dev/null" \
            % (vtyperPath, name, vtyperPath, name, primers)
        os.system(ePCR)
    # This populates vtyperResults with the verotoxin subtypes
    list1 = []
    if os.path.isfile("%s/%s.txt" % (vtyperPath, name)):
        ePCRresults = open("%s/%s.txt" % (vtyperPath, name), "r")
        for result in ePCRresults:
            # Only the lines without a # contain results
            if not "#" in result:
                count += 1
                # Split on \t
                data = result.split("\t")
                # The subtyping primer pair is the first entry on lines with results
                vttype = data[0].split("_")[0]
                # Push the name of the primer pair - stripped of anything after a _ to the dictionary
                if vttype not in list1:
                    list1.append(vttype)
    # Create a string of the entries in list1 joined with ";"
    string = ";".join(list1)
    metadata[name] = string
    # print json.dumps(metadata, sort_keys=True, indent=4, separators=(',', ': '))
    return metadata


def performBlast(name, path, targets, analysis, dictionary, refFilePath):
    """Runs the BLAST analysis and subsequent parsing of output files for GeneSeekr analysis"""
    output = {}
    geneSeekrList = []
    # GeneSeekr targets passed into the function
    for target in targets:
        # target is the full path and file name with extension - targetName removes path and extension
        targetName = os.path.split(target)[1].split(".")[0]
        # Set the out and query variables for the blastn command used below
        out = "%s/%s/geneSeekr/tmp/%s_%s.xml" % (path, name, name, targetName)
        query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
        # Perform a check to see if the files already exist, skip BLAST if they do
        if not os.path.isfile(out):
            # NcbiblastnCommandline runs a command line-like BLAST without actually invoking a system call
            blastn = NcbiblastnCommandline(query=query, db=target, outfmt=5, evalue=1e-10, out=out)
            stdout, stderr = blastn()
        # Parse the results
        result_handle = open(out)
        # Use the NCBIXML.parse functionality to rapidly parse the output files
        records = NCBIXML.parse(result_handle)
        # Counts the number of high-scoring pairs (HSPs) present in the output file
        numhsp = sum(line.count('<Hsp>') for line in iter(result_handle.readline, ""))
        # Skip if there are no HSPs
        if numhsp >= 1:
            # Since we scanned through result_handle looking for HSPs, the position of the read/write pointer
            # within the file is at the end. To reset it back to the beginning, .seek(0) is used
            result_handle.seek(0)
            # Go through each BLAST record
            for record in records:
                # print record.database_letters
                # Proceed iteratively through each alignment
                for alignment in record.alignments:
                    # And through each hsp in each alignment
                    for hsp in alignment.hsps:
                        # Since O-Typing requires a higher cutoff value in order to properly discriminate between
                        # o-types, a higher percent identity cutoff is specified
                        if "OType" in targetName:
                            percentIdentity = 1
                        # Otherwise, a cutoff of 80% should be fine
                        else:
                            percentIdentity = 0.8
                        # Invoke the percent identity cutoff. If the number of identities in the HSP is greater of equal
                        # to the total length of the subject sequence (multiplied by the cutoff), then the HSP passes
                        if hsp.identities >= alignment.length * percentIdentity:
                            # Calculate and format the percent identity of the match
                            PI = "%.2f" % (100 * (float(hsp.identities)/float(alignment.length)))
                            # Populate the dictionary with the results
                            # Note that since O-typing can have multiple 100% hits, I'm only taking one results - the
                            # way that the dictionary is populated, the H antigen is removed, so multiple results with
                            # different H antigens stop being unique
                            dictionary[name][targetName][alignment.title.split(" ")[0].split("_")[0]] = PI
                            if not alignment.title.split(" ")[0].split("_")[0] in geneSeekrList:
                                geneSeekrList.append(alignment.title.split(" ")[0].split("_")[0])
                            # Initiate the subtyping module
                            if "VT" in targetName:
                                output = vtyper(name, path, refFilePath)
                        # Otherwise, populate the dictionary indicating that the target was not found in the genome
                        else:
                            # Because O-typing has multiple targets per file, this populating of negative results would
                            # make a very bloated report - therefore don't populate negative results for each O-type
                            if not "OType" in targetName:
                                dictionary[name][targetName][alignment.title.split(" ")[0]] = "-"
        # If there are no HSPs, then there are no possible matches, so the target is not in the genome
        else:
            dictionary[name][targetName][targetName] = "-"
        # Close the file
        result_handle.close()
    # Create the report
    report = open("%s/%s/geneSeekr/%s_%sResults.tsv" % (path, name, name, analysis), "wb")
    # Populate appropriately
    report.write("Strain\t")
    # Get the targets into the header of the report
    for target in sorted(targets):
        report.write("%s\t" % os.path.split(target)[1].split(".")[0])
    # End the line
    report.write("\n")
    # Print strain name
    report.write("%s\t" % name)
    # Parse the dictionary to write presence/absence in the report
    for target in sorted(dictionary[name]):
        for tName in sorted(dictionary[name][target]):
            # print name, tName
            # Because the targetName of "OType" is "OType", the more useful tName (e.g. O100) is used instead
            # Additionally, the percent identity is included in brackets - will be 100% right now
            if target == "OType":
                report.write("%s (%s)\t" % (tName, dictionary[name][target][tName]))
            # tName is not required for non-OType genes, as the name of the gene is in the header, so only the
            # percent identity is written
            else:
                report.write("%s\t" % dictionary[name][target][tName])
    outputGeneSeekr = {}
    if geneSeekrList:
        string = ";".join(geneSeekrList)
    else:
        string = "N/A"
    outputGeneSeekr[name] = string
    # if output:
    return output, outputGeneSeekr

import threading
threadlock = threading.Lock()
# threadlock.acquire()
# threadlock.release()



def performMLST(name, path, refFilePath, genus):
    """Performs MLST analyses on strains with reference genomes that have a genus of Escherichia, Listeria
    or Salmonella"""
    metadata = {}
    # print json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
    # Set the percent identity required for a successful match to an allele
    percentIdentity = 0.98
    # As there are potentially multiple MLST schemes per organism, this checks to see
    # if multiple schemes exist
    numberSchemes = glob.glob("%s/Organism/%s/MLST/*" % (refFilePath, genus))
    # If there are not multiple schemes, then the alleles and profile folders should be present
    if not "alleles" in numberSchemes[0] and not "profile" in numberSchemes[0]:
        # Right now, this only gets the first MLST scheme - the ability to use multiple schemes can be added
        # TODO add the ability to use multiple schemes
        MLSTrefpath = numberSchemes[0]
    else:
        # Single schemed-organisms still need to populate the path
        MLSTrefpath = "%s/Organism/%s/MLST" % (refFilePath, genus)
    # Get the allele sequence file names
    alleles = glob.glob("%s/alleles/*.tfa" % MLSTrefpath)
    # Grab the profile file name
    profileFile = glob.glob("%s/profile/*.txt" % MLSTrefpath)
    # The gene names are present in the first line of the profile file
    # Note: if the profile files are ever updated, then the clonal complex column must be removed
    header = open(profileFile[0]).readline().rstrip().split("\t")
    # Get the MLST profiles for each sequence type
    with open(profileFile[0]) as profile:
        for line in profile:
            # MLSTcount will used to associate the gene name in header to the allele (e.g. adk 12)
            MLSTcount = 1
            # Don't need to get the header information again
            if not "ST" in line:
                # len(header) will be the same length as the data in line
                while MLSTcount < len(header):
                    # Remove newlines and split on tabs
                    data = line.rstrip().split("\t")
                    # Populate profileData with the sequence type, gene name, and the allele number
                    profileData[int(data[0])][header[MLSTcount]] = data[MLSTcount]
                    # Increment
                    MLSTcount += 1
    # Make the output path
    MLSToutpath = "%s/%s/MLST/tmp" % (path, name)
    make_path(MLSToutpath)
    # Perform the BLAST comparisons between all alleles of gene files, and the contigs file
    for allele in alleles:
        # Remove path and file extension information
        alleleName = os.path.split(allele)[1].split(".")[0]
        # Set the out and query variables prior to running BLAST
        out = "%s/%s_%s.xml" % (MLSToutpath, name, alleleName)
        query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
        # Perform a check to see if the files already exist, skip BLAST if they do
        if not os.path.isfile(out):
            # NcbiblastnCommandline runs a command line-like BLAST without actually invoking a system call
            blastn = NcbiblastnCommandline(query=query, db=allele, outfmt=5, evalue=1e-10, max_target_seqs=1, out=out)
            stdout, stderr = blastn()
        result_handle = open(out)
        # Use the NCBIXML.parse functionality to rapidly parse the output files
        records = NCBIXML.parse(result_handle)
        # Counts the number of high-scoring pairs (HSPs) present in the output file
        numhsp = sum(line.count('<Hsp>') for line in iter(result_handle.readline, ""))
        # Skip if there are no HSPs
        if numhsp >= 1:
            # Since we scanned through result_handle looking for HSPs, the position of the read/write pointer
            # within the file is at the end. To reset it back to the beginning, .seek(0) is used
            result_handle.seek(0)
            # Go through each BLAST record
            for record in records:
                # print record.database_letters
                # Proceed iteratively through each alignment
                for alignment in record.alignments:
                    # And through each hsp in each alignment
                    for hsp in alignment.hsps:
                        # Invoke the percent identity cutoff. If the number of identities in the HSP is greater of equal
                        # to the total length of the subject sequence (multiplied by the cutoff), then the HSP passes
                        if hsp.identities >= alignment.length * percentIdentity:
                            # Calculate and format the percent identity of the match
                            PI = "%.2f" % (100 * (float(hsp.identities)/float(alignment.length)))
                            # Populate the dictionary with the results
                            MLSTresults[name][alleleName][alignment.title.split(" ")[0].split("-")[1]] = PI
        else:
            MLSTresults[name][alleleName]["N/A"] = 0
        # Close the file
        result_handle.close()
    # print json.dumps(MLSTresults, sort_keys=True, indent=4, separators=(',', ': '))
    # Iterate through profileData and MLSTresults and find out how many matches to each sequence type occur
    # name is already defined in the function
    for alleleName in sorted(MLSTresults[name]):
        # Title is the allele number
        for title in MLSTresults[name][alleleName]:
            # Iterate through profileData now
            for sequenceType in profileData:
                # refallele is the allele number of the sequence type
                refAllele = profileData[sequenceType][alleleName]
                # If the allele number of the strain matches that of the sequence type...
                if title == refAllele:
                    # populate MLSTseqType with the appropriate information regarding that match
                    MLSTseqType[sequenceType][name][alleleName] = "%s (%s%%)" % (title, MLSTresults[name][alleleName][title])
    # print json.dumps(MLSTseqType, sort_keys=True, indent=4, separators=(',', ': '))
    # Iterate through MLSTseqType - if the number of matches between sequence type and the strain are equal to len(header)
    # (the total number of genes), then the strain is that sequence type
    # if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
    #     os.remove("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name))

    for sequenceType in sorted(MLSTseqType):
        for name in MLSTseqType[sequenceType]:
            # If there are seven matches out of a possible seven, then strain and sequence type are identical
            if len(MLSTseqType[sequenceType][name]) == len(header) - 1:
                # Write the match results to file
                threadlock.acquire()

                MLSTreport = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name), "wb")
                # Including name, as there will be a summary report where these data will be useful
                MLSTreport.write("Name\t")
                # Going back to header to get the gene names
                for entry in header:
                    # Expand ST to SequenceType for clarity
                    if "ST" in entry:
                        entry = "SequenceType"
                    # Write it to file
                    MLSTreport.write("%s\t" % entry)
                MLSTreport.write("\n")
                # Get the strain name, and sequence type information into the report
                MLSTreport.write("%s\t%s\t" % (name, sequenceType))
                # Populate the allele number information
                for alleleName in MLSTseqType[sequenceType][name]:
                    MLSTreport.write("%s\t" % MLSTseqType[sequenceType][name][alleleName])
                # MLSTreport.close()
                outputString = "(%s) %s" % (genus, sequenceType)
                metadata[name] = outputString
                # print name, outputString
                MLSTreport.close()
                threadlock.release()
                # metadata[name]["1.General"]["MLST_sequenceType"] = outputString
    fileSize = ""
    if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
        fileSize = os.stat("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name))
    if fileSize:
        if fileSize.st_size == 0:
            threadlock.acquire()
            # if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
            #     os.remove("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name))

            MLSTreport = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name), "wb")
            for name in MLSTresults:
                # print name
                metadata[name] = "N/A"
                # Write the match results to file
                # Including name, as there will be a summary report where these data will be useful
                MLSTreport.write("Name\t")
                # Going back to header to get the gene names
                for entry in header:
                    # Expand ST to SequenceType for clarity
                    if "ST" in entry:
                        entry = "SequenceType"
                    # Write it to file
                    MLSTreport.write("%s\t" % entry)
                MLSTreport.write("\n")
                # Get the strain name, and sequence type information into the report
                MLSTreport.write("%s\tN/A\t" % name)
                # Populate the allele number information
                for alleleName in MLSTresults[name]:
                    # print alleleName
                    for title, PI in MLSTresults[name][alleleName].iteritems():
                        MLSTreport.write("%s (%s%%)\t" % (title, PI))
            MLSTreport.close()
            threadlock.release()

    # Create a summary report of all MLST findings
    # Only attempt to append to the summary report if the strain has a MLST report
    if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
        # Because data from multiple files will be appended to the summary report, the header information
        # needs to be added first. If the file doesn't exist, then create it, and add header information
        if not os.path.isfile("%s/reports/%s_MLSTresults.tsv" % (path, genus)):
            # Same as above, so I will not add in-depth comments
            MLSTsummaryReport = open("%s/reports/%s_MLSTresults.tsv" % (path, genus), "wb")
            MLSTsummaryReport.write("Name\t")
            for entry in header:
                if "ST" in entry:
                    entry = "SequenceType"
                MLSTsummaryReport.write("%s\t" % entry)
            MLSTsummaryReport.write("\n")
            MLSTsummaryReport.close()
        # Get the results from the report - [1] means get the data from the second row
        MLSTsummary = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)).readlines()[1]
        metadata[name] = MLSTsummary.split("\t")[1]
        # Open the file again, but this time to append results
        MLSTsummaryReportData = open("%s/reports/%s_MLSTresults.tsv" % (path, genus), "ab")
        # Write results
        MLSTsummaryReportData.write("%s\n" % MLSTsummary)
        MLSTsummaryReportData.close()
    # # print json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
    if metadata:
        return metadata


def geneSeekrPrepProcesses(sampleName, path, assemblyMetadata, refFilePath):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction"""
    geneSeekrPrepArgs = []
    output = {}
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    if __name__ == 'geneSeekr':
        createGeneSeekrPool = Pool()
        # Prepare a tuple of the arguments (strainName and path)
        for name in sampleName:
            geneSeekrPrepArgs.append((name, path, assemblyMetadata, refFilePath))
        # This map function allows for multi-processing
        output = createGeneSeekrPool.map(runGeneSeekr, geneSeekrPrepArgs)
    return output


def runGeneSeekr((name, path, metadata, refFilePath)):
    """This function gathers relevant files, and calls appropriate functions - either GeneSeekr or MLST"""
    output = {}
    outputMLST = {}
    geneSeekrOutput = {}
    metadataOutput = defaultdict(make_dict)
    # Print a dot for each strain processed
    dotter()
    # All sequences are screened with the Illumina-specific UniVec database
    univecScreen(name, path, refFilePath)
    # Checks to ensure that the rMLST module has determined a closest reference genome - if not, then this
    # function will skip, as it is necessary to know the genus of the strain
    if os.path.isdir("%s/%s/referenceGenome" % (path, name)):
        # Parses the name of the reference genome to extract the genus
        genus = metadata[name]["1.General"]["referenceGenome"].split("_")[0]
        # Right now, only Escherichia, Listeria, and Salmonella genera can be run through the GeneSeekr/MLST functions
        if genus == "Escherichia" or genus == "Listeria" or genus == "Salmonella":
            # Grab the GeneSeekr
            targets = glob.glob("%s/Organism/%s/query_genes/*.fa" % (refFilePath, genus))
            # Make the directory
            make_path("%s/%s/geneSeekr/tmp" % (path, name))
            # Run GeneSeekr
            output, geneSeekrOutput = performBlast(name, path, targets, "geneSeekr", geneSeekrResults, refFilePath)
            # Run the MLST_typer
            outputMLST = performMLST(name, path, refFilePath, genus)
    # As there's no checks above to populate the dictionary if there are no results, the if statements below not only
    # populate the dictionary with the appropriate metadata headings, and return negative results
    if output:
        metadataOutput[name]["1.General"]["verotoxinProfile"] = output[name]
    else:
        metadataOutput[name]["1.General"]["verotoxinProfile"] = "N/A"

    if outputMLST:
        # print outputMLST[name]
        metadataOutput[name]["1.General"]["MLST_sequenceType"] = outputMLST[name]
    else:
        metadataOutput[name]["1.General"]["MLST_sequenceType"] = "N/A"

    if geneSeekrOutput:
        metadataOutput[name]["1.General"]["geneSeekrProfile"] = geneSeekrOutput[name]
    else:
        metadataOutput[name]["1.General"]["geneSeekrProfile"] = "N/A"
    return metadataOutput


def reportRemover(path):
    """As summary report files are created once, and then appended, if the pipeline is run multiple times, then the
    summary reports will contain results from each run. Reports should be cleared out each time the pipeline is run"""
    # Only these genera are processed with the MLST function
    genera = ["Escherichia", "Listeria", "Salmonella"]
    # Remove each report
    for genus in genera:
        if os.path.isfile("%s/reports/%s_MLSTresults.tsv" % (path, genus)):
            os.remove("%s/reports/%s_MLSTresults.tsv" % (path, genus))


def metadataFiller(metadata, geneSeekrList):
    """Properly populates the metadata dictionary - when I tried to populate the metadata dictionary within the
     multi-processed functions, it was returned as a list (likely due to how the files were returned. This essentially
     iterates through the list, and populates a dictionary appropriately"""
    # Make a copy of the metadata dictionary
    geneSeekrMetadata = metadata
    # Iterate through geneSeekrList
    for item in geneSeekrList:
        # each item in geneSeekrList is a dictionary entry
        # iterate through all the dictionaries
        for name in item:
            # The way the dictionaries were created, they should have the format:
            # metadataOutput[name]["1.General"]["geneSeekrProfile"] = geneSeekrOutput[name], so
            # e.g. "1.General"
            for generalCategory in item[name]:
                # e.g. "geneSeekrProfile"
                for specificCategory in item[name][generalCategory]:
                    # Populate the dictionary
                    geneSeekrMetadata[name][generalCategory][specificCategory] = str(item[name][generalCategory][specificCategory])
    # Return the beautifully-populated dictionary
    return geneSeekrMetadata


def functionsGoNOW(assembledFiles, path, assemblyMetadata, refFilePath):
    print "\nPerforming GeneSeekr analysis"
    # Clear out any summary reports from a previous iteration of the pipeline
    reportRemover(path)
    # Do everything - uniVec screening, geneSeeking, V-typing, and MLST analysis
    geneSeekrMetadataList = geneSeekrPrepProcesses(assembledFiles, path, assemblyMetadata, refFilePath)
    # print json.dumps(geneSeekrMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    geneSeekrMetadata = metadataFiller(assemblyMetadata, geneSeekrMetadataList)
    return geneSeekrMetadata