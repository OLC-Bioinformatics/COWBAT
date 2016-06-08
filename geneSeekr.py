import errno
import glob
import os
import shutil
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import jsonReportR
import metadataFiller
import threading

__author__ = 'adamkoziol'


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)


# Initialise the dictionary responsible for storing the report data
geneseekrresults = defaultdict(make_dict)
qualityresults = defaultdict(make_dict)
vtyperresults = defaultdict(make_dict)
profiledata = defaultdict(make_dict)
mlstresults = defaultdict(make_dict)
mlstseqtype = defaultdict(make_dict)

# Count variable is used in the dotter function
count = 0

# Initialise the threadlock
threadlock = threading.Lock()


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL.
    :param inpath: string of the supplied path
    """
    try:
        os.makedirs(inpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def univecscreen(name, path, reffilepath, command):
    """Runs all the sequences through the Illumina-specific UniVec database - output is a .tsv file
    in the path/name/univecscreen folder
    :param command:
    :param reffilepath:
    :param path:
    :param name: """
    univeccommand = defaultdict(make_dict)
    # Set the primer, query, and out variables for the BLAST search
    primers = "%s/illuminaPrimers/IlluminaPrimers.fa" % reffilepath
    out = "%s/%s/univecscreen/%s.tsv" % (path, name, name)
    query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
    # Make the path (if necessary)
    make_path("%s/%s/univecscreen" % (path, name))
    # If the analysis hasn't already been performed
    if not command[name]["UnivecCommand"]:
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
        blastn()
        univeccommand[name]["7.PipelineCommands"]["UnivecCommand"] = str(blastn)
    else:
        univeccommand[name]["7.PipelineCommands"]["UnivecCommand"] = command[name]["UnivecCommand"]
    return univeccommand


def vtyper(name, path, reffilepath, commands):
    """Uses epcr to subtype genomic verotoxin sequences
    :param commands:
    :param reffilepath:
    :param path:
    :param name:
    """
    metadata = defaultdict(make_dict)
    vtcommands = defaultdict(make_dict)
    # Because this is a multi-processed function, getting this dictionary to be properly formatted was difficult
    # in the end, I used if output, return output statements in the appropriate functions
    # global vtyperresults
    # Initialise count - this allows for the population of vtyperresults with unique values
    uniquecount = 0
    # Get the primers ready
    primers = "%s/vtyper/vtx_subtyping_primers.txt" % reffilepath
    # Make the output path
    vtyperpath = "%s/%s/vtyper" % (path, name)
    make_path(vtyperpath)
    # Copy the contigs file to the output path for ease of epcr indexing
    shutil.copyfile("%s/%s/%s_filteredAssembled.fasta" % (path, name, name),
                    "%s/%s/vtyper/%s_filteredAssembled.fasta" % (path, name, name))
    # For ease of processing set the contigs file as a variable
    query = "%s/%s_filteredAssembled.fasta" % (vtyperpath, name)
    # File checks to save re-performing processes
    if not commands[name]["v-typerFamapCommand"]:
        # famap and fahash are pre-processing steps performed on the query file
        famap = "famap -b %s/%s.famap %s 2>/dev/null" % (vtyperpath, name, query)
        os.system(famap)
        vtcommands[name]["7.PipelineCommands"]["v-typerFamapCommand"] = famap
    else:
        vtcommands[name]["7.PipelineCommands"]["v-typerFamapCommand"] = commands[name]["v-typerFamapCommand"]
    if not commands[name]["v-typerFahashCommand"]:
        fahash = "fahash -b %s/%s.hash %s/%s.famap 2>/dev/null" % (vtyperpath, name, vtyperpath, name)
        os.system(fahash)
        vtcommands[name]["7.PipelineCommands"]["v-typerFahashCommand"] = fahash
    else:
        vtcommands[name]["7.PipelineCommands"]["v-typerFahashCommand"] = commands[name]["v-typerFahashCommand"]
    if not commands[name]["v-typerEPCRCommand"]:
        # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
        # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
        # -m 10000 (Set variability for STS size for lookup), -n 1 (Set max allowed mismatches per primer for lookup)
        # -g 0 (Set max allowed indels per primer for lookup), -G (Print alignments in comments), -o {output file}
        epcr = "re-PCR -S %s/%s.hash -r + -m 10000 -n 1 -g 0 -G -q -o %s/%s.txt %s 2>/dev/null" \
               % (vtyperpath, name, vtyperpath, name, primers)
        os.system(epcr)
        vtcommands[name]["7.PipelineCommands"]["v-typerEPCRCommand"] = epcr
    else:
        vtcommands[name]["7.PipelineCommands"]["v-typerEPCRCommand"] = commands[name]["v-typerEPCRCommand"]
    # This populates vtyperresults with the verotoxin subtypes
    list1 = []
    if os.path.isfile("%s/%s.txt" % (vtyperpath, name)):
        epcrresults = open("%s/%s.txt" % (vtyperpath, name), "r")
        for result in epcrresults:
            # Only the lines without a # contain results
            if "#" not in result:
                uniquecount += 1
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
    # print json.dumps(vtcommands, sort_keys=True, indent=4, separators=(',', ': '))
    return metadata, vtcommands


def performblast(name, path, targets, analysis, dictionary, reffilepath, commands):
    """Runs the BLAST analysis and subsequent parsing of output files for GeneSeekr analysis
    :param commands:
    :param reffilepath:
    :param dictionary:
    :param analysis:
    :param targets:
    :param path:
    :param name:
    """
    output = {}
    vtcommands = defaultdict(make_dict)
    geneseekrcommands = defaultdict(make_dict)
    geneseekrlist = []
    if not commands[name]["GeneSeekrBlast"]:
        # GeneSeekr targets passed into the function
        for target in targets:
            # target is the full path and file name with extension - targetname removes path and extension
            targetname = os.path.split(target)[1].split(".")[0]
            # Set the out and query variables for the blastn command used below
            out = "%s/%s/geneSeekr/tmp/%s_%s.xml" % (path, name, name, targetname)
            query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
            # Perform a check to see if the files already exist, skip BLAST if they do
            # if not os.path.isfile(out):
            # NcbiblastnCommandline runs a command line-like BLAST without actually invoking a system call
            try:
                blastn = NcbiblastnCommandline(query=query, db=target, outfmt=5, evalue=1e-10, out=out)
                blastn()
                geneseekrcommands[name]["7.PipelineCommands"]["GeneSeekrBlast"] = str(blastn)
                # else:
                #     geneseekrcommands[name]["7.PipelineCommands"]["GeneSeekrBlast"] = commands[name]["GeneSeekrBlast"]
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
                                # Since O-Typing requires a higher cutoff value to properly discriminate between
                                # o-types, a higher percent identity cutoff is specified
                                if "OType" in targetname:
                                    percentidentity = 1
                                # Otherwise, a cutoff of 80% should be fine
                                else:
                                    percentidentity = 0.5
                                # Invoke the percent identity cutoff. If the number of identities in the HSP is greater 
                                # or equal to the total length of the subject sequence (multiplied by the cutoff), 
                                # then the HSP passes
                                if hsp.identities >= alignment.length * percentidentity:
                                    # Calculate and format the percent identity of the match
                                    pi = "%.2f" % (100 * (float(hsp.identities) / float(alignment.length)))
                                    # Populate the dictionary with the results
                                    # Note that since O-typing can have multiple 100% hits, I'm only taking one results 
                                    # the way that the dictionary is populated, the H antigen is removed, so multiple 
                                    # results with different H antigens stop being unique
                                    dictionary[name][targetname][alignment.title.split(" ")[0].split("_")[0]] = pi
                                    if not alignment.title.split(" ")[0].split("_")[0] in geneseekrlist:
                                        geneseekrlist.append(alignment.title.split(" ")[0].split("_")[0])
                                    # Initiate the subtyping module
                                    # In order to keep the program from initiating the subtyping module multiple times,
                                    # if the vtyper directory already exists for that strain, skip
                                    vtyperpath = "%s/%s/vtyper" % (path, name)
                                    if "VT" in targetname and not os.path.isdir(vtyperpath):
                                        output, vtcommands = vtyper(name, path, reffilepath, commands)
                                        # print name, path, reffilepath
                                # Otherwise, populate the dictionary indicating that the target 
                                # was not found in the genome
                                else:
                                    # Because O-typing has multiple targets per file, this populating of negative 
                                    # results would make a very bloated report - therefore don't populate negative 
                                    # results for each O-type
                                    if "OType" not in targetname:
                                        dictionary[name][targetname][alignment.title.split(" ")[0]] = "-"

                # If there are no HSPs, then there are no possible matches, so the target is not in the genome
                else:
                    dictionary[name][targetname][targetname] = "-"
                # Close the file
                result_handle.close()
            except TypeError:
                dictionary[name][targetname][targetname] = "-"
    else:
        geneseekrcommands[name]["7.PipelineCommands"]["GeneSeekrBlast"] = commands[name]["GeneSeekrBlast"]
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
            # Because the targetname of "OType" is "OType", the more useful tName (e.g. O100) is used instead
            # Additionally, the percent identity is included in brackets - will be 100% right now
            if target == "OType":
                report.write("%s (%s)\t" % (tName, dictionary[name][target][tName]))
            # tName is not required for non-OType genes, as the name of the gene is in the header, so only the
            # percent identity is written
            else:
                report.write("%s\t" % dictionary[name][target][tName])
    outputgeneseekr = {}
    if geneseekrlist:
        string = ";".join(geneseekrlist)
    else:
        string = "N/A"
    outputgeneseekr[name] = string
    # if output:
    return output, outputgeneseekr, vtcommands, geneseekrcommands


def performmlst(name, path, reffilepath, genus, commands):
    """Performs MLST analyses on strains with reference genomes that have a genus of Escherichia, Listeria
    or Salmonella
    :param commands:
    :param genus:
    :param reffilepath:
    :param path:
    :param name: """
    metadata = {}
    mlstcommands = defaultdict(make_dict)
    # print json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
    # Set the percent identity required for a successful match to an allele
    percentidentity = 0.98
    # As there are potentially multiple MLST schemes per organism, this checks to see
    # if multiple schemes exist
    numberschemes = glob.glob("%s/Organism/%s/MLST/*" % (reffilepath, genus))
    # If there are not multiple schemes, then the alleles and profile folders should be present
    if "alleles" not in numberschemes[0] and "profile" not in numberschemes[0]:
        # Right now, this only gets the first MLST scheme - the ability to use multiple schemes can be added
        # TODO add the ability to use multiple schemes
        mlstrefpath = numberschemes[0]
    else:
        # Single schemed-organisms still need to populate the path
        mlstrefpath = "%s/Organism/%s/MLST" % (reffilepath, genus)
    # Get the allele sequence file names
    alleles = glob.glob("%s/alleles/*.tfa" % mlstrefpath)
    # Grab the profile file name
    profilefile = glob.glob("%s/profile/*.txt" % mlstrefpath)
    # The gene names are present in the first line of the profile file
    # Note: if the profile files are ever updated, then the clonal complex column must be removed
    header = open(profilefile[0]).readline().rstrip().split("\t")
    # Get the MLST profiles for each sequence type
    with open(profilefile[0]) as profile:
        for line in profile:
            # mlstcount will used to associate the gene name in header to the allele (e.g. adk 12)
            mlstcount = 1
            # Don't need to get the header information again
            if "ST" not in line:
                # len(header) will be the same length as the data in line
                while mlstcount < len(header):
                    # Remove newlines and split on tabs
                    data = line.rstrip().split("\t")
                    # Populate profiledata with the sequence type, gene name, and the allele number
                    profiledata[int(data[0])][header[mlstcount]] = data[mlstcount]
                    # Increment
                    mlstcount += 1
    # Make the output path
    mlstoutpath = "%s/%s/MLST/tmp" % (path, name)
    make_path(mlstoutpath)
    # Perform the BLAST comparisons between all alleles of gene files, and the contigs file
    for allele in alleles:
        # Remove path and file extension information
        allelename = os.path.split(allele)[1].split(".")[0]
        # Set the out and query variables prior to running BLAST
        out = "%s/%s_%s.xml" % (mlstoutpath, name, allelename)
        query = "%s/%s/%s_filteredAssembled.fasta" % (path, name, name)
        # Perform a check to see if the files already exist, skip BLAST if they do
        if not os.path.isfile(out):
            # NcbiblastnCommandline runs a command line-like BLAST without actually invoking a system call
            blastn = NcbiblastnCommandline(query=query, db=allele, outfmt=5, evalue=1e-10, max_target_seqs=1, out=out)
            blastn()
            mlstcommands[name]["7.PipelineCommands"]["MLSTBlastCommand"] = str(blastn)
        else:
            mlstcommands[name]["7.PipelineCommands"]["MLSTBlastCommand"] = commands[name]["MLSTBlastCommand"]
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
                        if hsp.identities >= alignment.length * percentidentity:
                            # Calculate and format the percent identity of the match
                            pi = "%.2f" % (100 * (float(hsp.identities) / float(alignment.length)))
                            # Populate the dictionary with the results
                            mlstresults[name][allelename][alignment.title.split(" ")[0].split("-")[1]] = pi
        else:
            mlstresults[name][allelename]["N/A"] = 0
        # Close the file
        result_handle.close()
    # print json.dumps(mlstresults, sort_keys=True, indent=4, separators=(',', ': '))
    # Iterate through profiledata and mlstresults and find out how many matches to each sequence type occur
    # name is already defined in the function
    for allelename in sorted(mlstresults[name]):
        # Title is the allele number
        for title in mlstresults[name][allelename]:
            # Iterate through profiledata now
            for sequencetype in profiledata:
                # refallele is the allele number of the sequence type
                refallele = profiledata[sequencetype][allelename]
                # If the allele number of the strain matches that of the sequence type...
                if title == refallele:
                    # populate mlstseqtype with the appropriate information regarding that match
                    mlstseqtype[sequencetype][name][allelename] = "%s (%s%%)" % (
                        title, mlstresults[name][allelename][title])
    # print json.dumps(mlstseqtype, sort_keys=True, indent=4, separators=(',', ': '))
    # Iterate through mlstseqtype. If the number of matches between sequence type and strain are equal to len(header)
    # (the total number of genes), then the strain is that sequence type

    for sequencetype in sorted(mlstseqtype):
        for name in mlstseqtype[sequencetype]:
            # If there are seven matches out of a possible seven, then strain and sequence type are identical
            if len(mlstseqtype[sequencetype][name]) == len(header) - 1:
                # Write the match results to file
                threadlock.acquire()

                mlstreport = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name), "wb")
                # Including name, as there will be a summary report where these data will be useful
                mlstreport.write("Name\t")
                # Going back to header to get the gene names
                for entry in header:
                    # Expand ST to SequenceType for clarity
                    if "ST" in entry:
                        entry = "SequenceType"
                    # Write it to file
                    mlstreport.write("%s\t" % entry)
                mlstreport.write("\n")
                # Get the strain name, and sequence type information into the report
                mlstreport.write("%s\t%s\t" % (name, sequencetype))
                # Populate the allele number information
                for allelename in sorted(mlstseqtype[sequencetype][name]):
                    # print allelename, mlstseqtype[sequencetype][name][allelename]
                    mlstreport.write("%s\t" % mlstseqtype[sequencetype][name][allelename])
                # mlstreport.close()
                outputstring = "(%s) %s" % (genus, sequencetype)
                metadata[name] = outputstring
                # print name, outputstring
                mlstreport.close()
                threadlock.release()
                # print name, outputstring
    # for sequencetype in sorted(mlstseqtype):
    if name not in metadata:
        threadlock.acquire()
        mlstreport = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name), "wb")
        mlstreport.write("Name\t")
        for entry in header:
            # Expand ST to SequenceType for clarity
            if "ST" in entry:
                entry = "SequenceType"
                # Write it to file
            mlstreport.write("%s\t" % entry)
        mlstreport.write("\n")
        mlstreport.write("%s\tN/A\t" % name)
        for allelename in sorted(mlstresults[name]):
            # Going back to header to get the gene names
            # Title is the allele number
            for title in sorted(mlstresults[name][allelename]):
                if title:
                    mlstreport.write("%s (%s%%)\t" % (title, mlstresults[name][allelename][title]))
        mlstreport.close()
        threadlock.release()
        # print name, len(mlstseqtype[sequencetype][name])
        # metadata[name]["1.General"]["MLST_sequenceType"] = outputstring
    # print name
    filesize = ""
    if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
        filesize = os.stat("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name))
    if filesize:
        if filesize.st_size == 0:
            threadlock.acquire()
            # if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
            #     os.remove("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name))
            mlstreport = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name), "wb")
            for name in mlstresults:
                print name
                metadata[name] = "N/A"
                # Write the match results to file
                # Including name, as there will be a summary report where these data will be useful
                mlstreport.write("Name\t")
                # Going back to header to get the gene names
                for entry in header:
                    # Expand ST to SequenceType for clarity
                    if "ST" in entry:
                        entry = "SequenceType"
                    # Write it to file
                    mlstreport.write("%s\t" % entry)
                mlstreport.write("\n")
                # Get the strain name, and sequence type information into the report
                mlstreport.write("%s\tN/A\t" % name)
                # Populate the allele number information
                for allelename in mlstresults[name]:
                    # print allelename
                    for title, pi in mlstresults[name][allelename].iteritems():
                        mlstreport.write("%s (%s%%)\t" % (title, pi))
            mlstreport.close()
            threadlock.release()

    # Create a summary report of all MLST findings
    # Only attempt to append to the summary report if the strain has a MLST report
    if os.path.isfile("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)):
        # Because data from multiple files will be appended to the summary report, the header information
        # needs to be added first. If the file doesn't exist, then create it, and add header information
        if not os.path.isfile("%s/reports/%s_MLSTresults.tsv" % (path, genus)):
            # Same as above, so I will not add in-depth comments
            mlstsummaryreport = open("%s/reports/%s_MLSTresults.tsv" % (path, genus), "wb")
            mlstsummaryreport.write("Name\t")
            for entry in header:
                if "ST" in entry:
                    entry = "SequenceType"
                mlstsummaryreport.write("%s\t" % entry)
            mlstsummaryreport.write("\n")
            mlstsummaryreport.close()
        # Get the results from the report - [1] means get the data from the second row
        mlstsummary = open("%s/%s/MLST/%s_MLSTreport.tsv" % (path, name, name)).readlines()[1]
        metadata[name] = mlstsummary.split("\t")[1]
        # Open the file again, but this time to append results
        mlstsummaryreportdata = open("%s/reports/%s_MLSTresults.tsv" % (path, genus), "ab")
        # Write results
        mlstsummaryreportdata.write("%s\n" % mlstsummary)
        mlstsummaryreportdata.close()
    # # print json.dumps(metadata[name], sort_keys=True, indent=4, separators=(',', ': '))
    if metadata:
        return metadata, mlstcommands
    else:
        return mlstcommands


def geneseekrprepprocesses(samplename, path, assemblymetadata, reffilepath, commands):
    """A helper function to make a pool of processes to allow for a multi-processed approach to error correction
    :param commands:
    :param reffilepath:
    :param assemblymetadata:
    :param path:
    :param samplename:
    """
    geneseekrprepargs = []
    # This used to check to see if the __name__ == '__main__', but since this is run as a module, the __name__
    # must equal the name of the script
    creategeneseekrpool = Pool()
    # Prepare a tuple of the arguments (strainName and path)
    for name in samplename:
        geneseekrprepargs.append((name, path, assemblymetadata, reffilepath, commands))
    # This map function allows for multi-processing
    output = creategeneseekrpool.map(rungeneseekr, geneseekrprepargs)
    return output


def rungeneseekr((name, path, metadata, reffilepath, commands)):
    """This function gathers relevant files, and calls appropriate functions - either GeneSeekr or MLST"""
    output = {}
    outputmlst = {}
    geneseekroutput = {}
    vtcommands = {}
    mlstcommands = {}
    geneseekrcommands = {}
    metadataoutput = defaultdict(make_dict)
    # Print a dot for each strain processed
    dotter()
    # All sequences are screened with the Illumina-specific UniVec database
    univeccommand = univecscreen(name, path, reffilepath, commands)
    # Checks to ensure that the rMLST module has determined a closest reference genome - if not, then this
    # function will skip, as it is necessary to know the genus of the strain
    if os.path.isdir("%s/%s/referenceGenome" % (path, name)):
        # Parses the name of the reference genome to extract the genus
        genus = metadata[name]["1.General"]["referenceGenome"].split("_")[0]
        # Right now, only Escherichia, Listeria, and Salmonella genera can be run through the GeneSeekr/MLST functions
        if genus == "Escherichia" or genus == "Listeria" or genus == "Salmonella":
            # Grab the GeneSeekr
            targets = glob.glob("%s/Organism/%s/query_genes/*.fa" % (reffilepath, genus))
            # Make the directory
            make_path("%s/%s/geneSeekr/tmp" % (path, name))
            # Run GeneSeekr
            output, geneseekroutput, vtcommands, geneseekrcommands = performblast(name, path, targets, "geneSeekr",
                                                                                  geneseekrresults, reffilepath,
                                                                                  commands)
            # Run the MLST_typer
            try:
                outputmlst, mlstcommands = performmlst(name, path, reffilepath, genus, commands)
            except KeyError:
                pass
    # print name, geneseekroutput, vtcommands, geneseekrcommands, outputmlst, mlstcommands
    # As there's no checks above to populate the dictionary if there are no results, the if statements below not only
    # populate the dictionary with the appropriate metadata headings, and return negative results
    try:
        if output:
            metadataoutput[name]["1.General"]["verotoxinProfile"] = output[name]
        else:
            metadataoutput[name]["1.General"]["verotoxinProfile"] = "N/A"
    except KeyError:
        print name, "verotoxinProfile"
    try:
        if outputmlst:
            # print outputmlst[name]
            metadataoutput[name]["1.General"]["MLST_sequenceType"] = outputmlst[name]
        else:
            metadataoutput[name]["1.General"]["MLST_sequenceType"] = "N/A"
    except KeyError:
        metadataoutput[name]["1.General"]["MLST_sequenceType"] = "N/A"
    try:
        if geneseekroutput:
            metadataoutput[name]["1.General"]["geneSeekrProfile"] = geneseekroutput[name]
        else:
            metadataoutput[name]["1.General"]["geneSeekrProfile"] = "N/A"
    except KeyError:
        print name, "geneSeekrProfile"
    try:
        if vtcommands:
            metadataoutput[name]["7.PipelineCommands"]["v-typerEPCRCommand"] = vtcommands[name]["7.PipelineCommands"][
                "v-typerEPCRCommand"]
            metadataoutput[name]["7.PipelineCommands"]["v-typerFahashCommand"] = vtcommands[name]["7.PipelineCommands"][
                "v-typerFahashCommand"]
            metadataoutput[name]["7.PipelineCommands"]["v-typerFamapCommand"] = vtcommands[name]["7.PipelineCommands"][
                "v-typerFamapCommand"]
        else:
            metadataoutput[name]["7.PipelineCommands"]["v-typerEPCRCommand"] = commands[name]["v-typerEPCRCommand"]
            metadataoutput[name]["7.PipelineCommands"]["v-typerFahashCommand"] = commands[name]["v-typerFahashCommand"]
            metadataoutput[name]["7.PipelineCommands"]["v-typerFamapCommand"] = commands[name]["v-typerFamapCommand"]
    except KeyError:
        print name, "typerEPCRCommand"
    try:
        if geneseekrcommands:
            metadataoutput[name]["7.PipelineCommands"]["GeneSeekrBlast"] = geneseekrcommands[name]["7.PipelineCommands"][
                "GeneSeekrBlast"]
        else:
            metadataoutput[name]["7.PipelineCommands"]["GeneSeekrBlast"] = commands[name]["GeneSeekrBlast"]
    except KeyError:
        print name, "GeneSeekrBlast"
    try:
        if univeccommand:
            metadataoutput[name]["7.PipelineCommands"]["UnivecCommand"] = univeccommand[name]["7.PipelineCommands"][
                "UnivecCommand"]
        else:
            metadataoutput[name]["7.PipelineCommands"]["UnivecCommand"] = commands[name]["UnivecCommand"]
    except KeyError:
        print name, "UnivecCommand"
    try:
        if mlstcommands:
            metadataoutput[name]["7.PipelineCommands"]["MLSTBlastCommand"] = mlstcommands[name]["7.PipelineCommands"][
                "MLSTBlastCommand"]
        else:
            metadataoutput[name]["7.PipelineCommands"]["MLSTBlastCommand"] = commands[name]["MLSTBlastCommand"]
    except KeyError:
        print name, "MLSTBlastCommand"

    return metadataoutput


def reportremover(path):
    """As summary report files are created once, and then appended, if the pipeline is run multiple times, then the
    summary reports will contain results from each run. Reports should be cleared out each time the pipeline is run
    :param path: """
    # Only these genera are processed with the MLST function
    genera = ["Escherichia", "Listeria", "Salmonella"]
    # Remove each report
    for genus in genera:
        if os.path.isfile("%s/reports/%s_MLSTresults.tsv" % (path, genus)):
            os.remove("%s/reports/%s_MLSTresults.tsv" % (path, genus))


def functionsgonow(assembledfiles, path, assemblymetadata, reffilepath, commands):
    print "\nPerforming GeneSeekr analysis"
    # Clear out any summary reports from a previous iteration of the pipeline
    reportremover(path)
    # Do everything - uniVec screening, geneSeeking, V-typing, and MLST analysis
    geneseekrmetadatalist = geneseekrprepprocesses(assembledfiles, path, assemblymetadata, reffilepath, commands)
    # print json.dumps(geneseekrmetadata, sort_keys=True, indent=4, separators=(',', ': '))
    geneseekrmetadata = metadataFiller.filler(assemblymetadata, geneseekrmetadatalist)
    jsonReportR.jsonR(assembledfiles, path, geneseekrmetadata, "Collection")
    return geneseekrmetadata
