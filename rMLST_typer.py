__author__ = 'akoziol'

""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
import subprocess, os, glob, time, sys, shlex, re, threading, json, mmap, shutil, errno
import jsonReportR
# Initialise variables
count = 0
dqueue = Queue()
blastqueue = Queue()
parsequeue = Queue()
testqueue = Queue()
plusqueue = Queue()
plusdict = {}
genedict = defaultdict(list)
blastpath = {}
threadlock = threading.Lock()
genomes = []
scheme = []
header = []
body = []


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
sequenceTypes = defaultdict(make_dict)
strainTypes = defaultdict(make_dict)
mismatch = defaultdict(make_dict)

sequenceTypesMLST = defaultdict(make_dict)
strainTypesMLST = defaultdict(make_dict)
mismatchMLST = defaultdict(make_dict)
genomeCovered = defaultdict(make_dict)
testDict = defaultdict(make_dict)

profileMismatch = defaultdict(make_dict)

def dotter():
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
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


def makeblastdb(dqueue):
    while True:  # while daemon
        fastapath = dqueue.get() # grabs fastapath from dqueue
        nhr = "%s.nhr" % (fastapath)
        if not os.path.isfile(str(nhr)):
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, fastapath)))
            dotter()
        dqueue.task_done() # signals to dqueue job is done
        sys.exit()


def makedbthreads(fastas):
    ''' Setup and create threads for class'''
    for i in range(len(fastas)):
        threads = Thread(target=makeblastdb, args=(dqueue,))
        threads.setDaemon(True)
        threads.start()
    for fasta in fastas:
        dqueue.put(fasta)
    dqueue.join() #wait on the dqueue until everything has been processed


def xmlout (fasta, genome):
    """Extracts the name of the file by stripping off extensions"""
    gene = re.search('\/(\w+)\.tfa', fasta)
    path = re.search('(.+)\/(.+)\/(.+?)\.fa', genome)
    return path, gene


class runblast(threading.Thread):
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)

    def run(self):
        while True:
            global blastpath, plusdict, genomeCovered
            genome, fasta, blastexist = self.blastqueue.get()
            path, gene = xmlout(fasta, genome)
            out = "%s/tmp/%s.%s.xml" % (path.group(1), path.group(3), gene.group(1))
            threadlock.acquire()
            blastpath[out] = {path.group(3): (gene.group(1),)}
            threadlock.release()
            if not os.path.isfile(out):
                dotter()
                blastn = NcbiblastnCommandline(query=genome, db=fasta, evalue=1e-40, out=out, outfmt=5, perc_identity=100)
                stdout, stderr = blastn()
            if not any(blastpath):
                print out
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes):
    '''Setup and create  threads for blastn and xml path'''
    blastexist = {}
    for i in range(len(fastas)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        threads.start()
    for genome in genomes:
        for fasta in fastas:
            blastqueue.put((genome, fasta, blastexist))
        blastqueue.join()


class blastparser(threading.Thread):  # records, genomes):
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)

    def rMLSTadder(self):
        xml, genomes, mm, num = self.parsequeue.get()
        print xml, genomes, mm, num

    def run(self):
        while True:
            global plusdict, genedict
            col1 = []
            alignmentLength = 0
            alleleMatches = 0
            xml, genomes, mm, num = self.parsequeue.get()
            records = NCBIXML.parse(mm)
            numhsp = sum(line.count('<Hsp>') for line in iter(mm.readline, ""))
            if numhsp >= 1:
                mm.seek(0)
                col = 'N'
                for record in records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            threadlock.acquire()  # precaution
                            if hsp.identities == alignment.length:
                                col = alignment.title.split('_')[-1]  # MLST type
                                alignmentLength += alignment.length
                                alleleMatches += 1
                                col = re.sub(" No definition line", "", col)
                                col1.append(col)
                                # print alignment.title.split(" ")[0]
                            threadlock.release()  # precaution for populate dictionary with GIL
                # print len(col1)
                # if len(col1) > 1:
                # if col1:
                    # print " ".join(col1)
                col1.sort()
                mm.seek(0)
                records = NCBIXML.parse(mm)
                if not col1:
                    # print col
                    for record in records:
                        # print record
                        for alignment in record.alignments:
                            for hsp in alignment.hsps:
                                # threadlock.acquire()  # precaution
                                # if col == 'N' and alignment.length > hsp.identities >= alignment.length * 0.98:
                                    # col = 'B'
                                    # print mm, num, alignment.title.split(" ")[0]
                                    # if hsp.identities >= alignment.length * 0.3:
                                    #     print genomes, hsp.identities, alignment.length
                                        # PI = hsp.identities / alignment.length
                                        # print genomes, num, "%.2f" % (PI)
                                        # self.rMLSTadder()
                                for genome in genomes:
                                    for gene in genomes[genome]:
                                        # print gene, genome, alignment.title.split(" ")[0], hsp.identities, alignment.length
                                        if alignment.length > hsp.identities >= alignment.length * 0.95:
                                        # if gene == "BACT000048":
                                        #     PI = (hsp.identities / alignment.length) * 100
                                        #     print gene, genome, alignment.title.split(" ")[0], hsp.identities, alignment.length
                                            colPI = (hsp.identities / alignment.length) * 100
                                            if genome not in plusdict:
                                                plusdict[genome] = defaultdict(str)
                                            if gene[-2:] not in plusdict[genome]:
                                                plusdict[genome][gene[-2:]] = colPI
                for genome in genomes:
                    # if alignmentLength > 0:
                        # print genome, alignmentLength
                    # print genome
                    name = re.split("_filteredAssembled", str(genome))[0]
                    # for aL in genomeCovered[name]:
                        # aL += alignmentLength

                        # genomeCovered[str(genome)] = aL
                    for gene in genomes[genome]:
                        # testDict[name][gene][col] = alignmentLength
                        # genomeCovered[name] += alignmentLength
                        # if alignment.length > hsp.identities >= alignment.length * 0.95:
                        # if gene == "BACT000048":
                        #     PI = (hsp.identities / alignment.length) * 100
                        #     print gene, genome, alignment.title.split(" ")[0], hsp.identities, alignment.length
                        if genome not in plusdict:
                            plusdict[genome] = defaultdict(str)
                        if gene[-2:] not in plusdict[genome]:
                            if len(col1) > 1:
                                plusdict[genome][gene[-2:]] = " ".join(col1)
                                # testDict[name][gene][" ".join(col1)][alleleMatches] = alignmentLength
                                testDict[name] += alleleMatches
                                genomeCovered[name] += alignmentLength
                            else:
                                plusdict[genome][gene[-2:]] = col
                                # testDict[name][gene][" ".join(col1)][alleleMatches] = alignmentLength
                                testDict[name] += alleleMatches
                                genomeCovered[name] += alignmentLength
                                # threadlock.release()  # precaution for populate dictionary with GIL
            else:
                for genome in genomes:
                    for gene in genomes[genome]:
                        if genome not in plusdict:
                            plusdict[genome] = defaultdict(str)
                        if gene[-2:] not in plusdict[genome]:
                            plusdict[genome][gene[-2:]] = "N"
            # print alignmentLength
            # print json.dumps(genomeCovered, sort_keys=True, indent=4, separators=(',', ': '))
            dotter()
            mm.close()
            self.parsequeue.task_done()
            # return genomeCovered




def parsethreader(bPath, Genomes):
    """Sets up multithreaded blast"""
    global plusdict
    # global genomeCovered
    dotter()
    for i in range(len(Genomes)):
        threads = blastparser(parsequeue)
        threads.setDaemon(True)
        threads.start()
    progress = len(bPath)
    for xml in bPath:
        handle = open(xml, 'r')
        mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        handle.close()
        parsequeue.put((xml, bPath[xml], mm, progress))
        parsequeue.join()


def blaster(markers, strains, path, out, experimentName, refFilesPath):
    '''
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    ALL PATHS REQUIRE TRAILING SLASHES!!!
    '''
    global count, genedict, blastpath
    #retrieve rMLST markers from input
    fastas = glob.glob(markers + "*.tfa")
    # #retrieve genomes from input
    for name in strains:
        genomeFile = glob.glob("%s/%s/*_filteredAssembled.fasta" % (path, name))
        # print genomeFile
        if genomeFile:
            genomes.append(genomeFile[0])
        # else:
        #     pass
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    if os.path.isfile('%s/blastxmldict.json' % path):
        blastpath = json.load(open('%s/blastxmldict.json' % path))
    else:
        # make blastn threads and retrieve xml file locations
        blastnthreads(fastas, genomes)
        # , indent=4, separators=(',', ': ')
        json.dump(blastpath, open('%s/blastxmldict.json' % path, 'wb'), sort_keys=True, indent=4)
    print "\n[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
    parsethreader(blastpath, fastas)
    csvheader = 'Strain'
    row = ""
    rowcount = 0
    for genomerow in plusdict:
        row += "\n" + genomerow
        rowcount += 1
        for generow in sorted(plusdict[genomerow]):
            if rowcount <= 1:
                csvheader += ', BACT0000' + generow.zfill(2)
            row += ',' + str(plusdict[genomerow][generow])
    make_path("%s/reports" % path)
    with open("%s/reports/%s_rMLSTresults.tsv" % (out, experimentName), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    return plusdict


def determineReferenceGenome(plusdict, path, metadata, refFilesPath):
    """Find the reference genome with the closest number of identical alleles to each
    sequenced strain"""
    # Retrieve the rMLST profiles of the reference genomes
    with open("%s/referenceGenomes/referenceGenome_rMLSTprofiles.bak" % refFilesPath, "r") as referenceFile:
        sequenceTypes = json.load(referenceFile)
        referenceFile.close()
    print "\nDetermining closest reference genome."
    for genome in plusdict:
        # if not os.path.isdir("%s/%s")
        dotter()
        bestCount = 0
        for sType in sorted(sequenceTypes):
            count = 0
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
                # print genome, sType, bactNum, allele
                # bactNum is only the two trailing digits of the gene name
                fullBact = "BACT0000%s" % bactNum
                # Checks to see if the key (e.g. BACT000001) is in the dictionary of reference genomes
                if str(fullBact) in sequenceTypes[sType]:
                    # Assigns the corresponding allele to refAllele
                    referenceAllele = sequenceTypes[sType][fullBact]
                    # Part of the reference genome rMLST profile includes genes that match two alleles
                    # if this is encountered, split on the space, and check each allele separately
                    if re.search(" ", referenceAllele):
                        splitAllele = referenceAllele.split(" ")
                        # If the alleles match, increment the count
                        if allele == splitAllele[0] or allele == splitAllele[1]:
                            count += 1
                        # If not, add the mismatch details to the mismatch dictionary
                        else:
                            mismatch[fullBact][allele] = sequenceTypes[sType][fullBact]
                    # No spaces encountered in the allele definitions
                    else:
                        if allele == referenceAllele:
                            count += 1
                        else:
                            mismatch[fullBact][allele] = sequenceTypes[sType][fullBact]
                # If the key isn't in the dictionary, add a "N" to the mismatch dictionary
                else:
                    mismatch[fullBact][allele] = "N"
            # In order to find the best match, compare the score of the previous best match
            # to this current match, if this new match is better, then proceed
            if count > bestCount:
                # Set the best count to current count
                bestCount = count
                # Clear out the previous information from the dictionary
                strainTypes[genome].clear()
                # If mismatches are recorded in the dictionary
                if mismatch:
                    # Traverse the dictionary
                    for geneName in sorted(mismatch):
                        for observedAllele, refAllele in sorted(mismatch[geneName].iteritems()):
                            # Populate the metadata dictionary
                            strainTypes[genome][sType]["gene"][geneName]["observedAllele"][observedAllele]["referenceAllele"] = refAllele
                            strainTypes[genome][sType]["NumIdenticalAlleles"] = bestCount
                # If no mismatches, populate the dictionary with the strain name, reference genome name, and the count
                else:
                    strainTypes[genome][sType]["NumIdenticalAlleles"] = bestCount
            # Clear the dictionary for the next round of analysis
            mismatch.clear()
    # Populate the metadata dictionary with the appropriate data
    for strain in sorted(strainTypes):
        strainTrimmed = re.split("_filteredAssembled", strain)[0]
        make_path("%s/%s/referenceGenome" % (path, strainTrimmed))
        for reference in sorted(strainTypes[strain]):
            # Copy the reference genome to the referenceGenome subfolder in the strain directory
            shutil.copyfile("%s/referenceGenomes/%s.fasta" % (refFilesPath, reference), "%s/%s/referenceGenome/%s.fasta" % (path, strainTrimmed, reference))

            metadata[strainTrimmed]["6.rMLSTmatchestoRef"] = strainTypes[strain][reference]
            # print strainTypes[strain][reference]
            # metadata[strainTrimmed]["6.rMLSTmatchestoRef"]["NumIdenticalAlleles"] = bestCount
    return metadata


def profileUpdater(plusdict, path, metadata, refFilesPath, genomes, lastEntry, strainTypesMLST):
    print "\nAdding rMLST profiles."
    testDictionary = {}
    # new_dict = defaultdict(list)
    new_list = []
    listr = []
    for genome in plusdict:
        # print genome, plusdict[genome]
        # if plusdict[genome] not in testDictionary:
        #     print genome
            # testDictionary[genome]
        # print genome
        for bactNum, allele in sorted(plusdict[genome].iteritems()):
        #     testDictionary
        #     comparison = list(set(allele))
    # print comparison
            listr.append([str(bactNum), str(allele)])
        testDictionary[genome] = listr
        #
        # print genome, listr
        listr = []
    # for genome in testDictionary:
    #     print genome, testDictionary[genome]
        # print json.dumps(testDictionary, sort_keys=True, indent=4, separators=(',', ': '))
    for k, v in testDictionary.iteritems():
    #     print k, v
        if v not in new_list:
            # print k, v
            new_list.append(v)
        # else:
        #     print k
    profile = open("%s/rMLST/profile/rMLST_scheme.txt" % refFilesPath, "ab")

    for i in new_list:
        lastEntry += 1
        if lastEntry < 1000000:
            lastEntry = 1000000
        profile.write("%s_CFIA\t" % lastEntry)
        for j in i:
            # print j
            blarg = j[1]
            # print lastEntry, j, blarg
            profile.write("%s\t" % blarg)
        profile.write("\n")
    profile.close()

    profile = open("%s/rMLST/profile/rMLST_scheme.txt" % refFilesPath)
    header = []
    body = []
    for line in profile:
        # Skip the first line of the scheme
        if re.search("rST", line):
            subline = line.split("\t")
            for subsubline in subline:
                if re.search("BACT", subsubline):
                    header.append(subsubline)
        else:
            subline = line.split("\t")
            body.append(subline)
    profile.close()

    sequenceTypesMLST.clear()
    #
    for line in body:
        for i in range(len(header)):
            # print line[0], header[i], line[i + 1]
            sequenceTypesMLST[line[0]][header[i]] = line[i + 1]
    for genome in plusdict:
        dotter()
        bestCount = 0
        for sType in sequenceTypesMLST:
            count = 0
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
                fullBact = "BACT0000%s" % bactNum
                if allele == sequenceTypesMLST[sType][fullBact]:
                    count += 1
            if count > bestCount:
                bestCount = count
                strainTypesMLST[genome].clear()
                strainTypesMLST[genome][str(sType)]["rMLSTIdenticalAlleles"] = bestCount
    return strainTypesMLST


def determineSubtype(plusdict, path, metadata, refFilesPath):
    """Same as above (determineReferenceGenome), but this determines the rMLST sequence type
    Comments are as above, as this is very similar"""
    profile = open("%s/rMLST/profile/rMLST_scheme.txt" % refFilesPath)
    lastEntry = int(open("%s/rMLST/profile/rMLST_scheme.txt" % refFilesPath).readlines()[-1].split("\t")[0].split("_CFIA")[0])
    # print lastEntry
    for line in profile:
        # Skip the first line of the scheme
        if re.search("rST", line):
            subline = line.split("\t")
            for subsubline in subline:
                if re.search("BACT", subsubline):
                    header.append(subsubline)
        else:
            subline = line.split("\t")
            body.append(subline)
    profile.close()
    for line in body:
        for i in range(len(header)):
            # print line[0], header[i], line[i + 1]
            sequenceTypesMLST[line[0]][header[i]] = line[i + 1]
    print "\nDetermining sequence types."
    for genome in plusdict:
        dotter()
        bestCount = 0
        for sType in sequenceTypesMLST:
            count = 0
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
                fullBact = "BACT0000%s" % bactNum
                if allele == sequenceTypesMLST[sType][fullBact]:
                    count += 1
            #
            if count > bestCount:
                bestCount = count
                strainTypesMLST[genome].clear()
                strainTypesMLST[genome][str(sType)]["rMLSTIdenticalAlleles"] = bestCount
        #
        if bestCount < len(header):
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
                # fullBact = "BACT0000%s" % bactNum
                profileMismatch[genome][bactNum] = allele
    # Somehow, there are some files which do not have a single allele in common with any
    # of the profiles present in the rMLST scheme. These are addressed here.
    strainTypesMLSTprofiled = profileUpdater(profileMismatch, path, metadata, refFilesPath, genomes, lastEntry, strainTypesMLST)
    # print json.dumps(strainTypesMLSTprofiled, sort_keys=True, indent=4, separators=(',', ': '))
    for genome in plusdict:
        # If the genome key is not in the strainTypesMLST dictionary
        if genome not in sorted(strainTypesMLSTprofiled):
            strainTrimmed = re.split("_filteredAssembled", genome)[0]
            metadata[strainTrimmed]["5.rMLST"]["rMLSTSequenceType"] = "N/A"
    # Populates the metadata dictionary with data from strains present in the strainTypesMLST dictionary
    for strain in sorted(strainTypesMLSTprofiled):
        strainTrimmed = re.split("_filteredAssembled", strain)[0]
        for INTreference in sorted(strainTypesMLSTprofiled[strain]):
            reference = str(INTreference)
            metadata[strainTrimmed]["5.rMLST"] = strainTypesMLSTprofiled[strain][reference]
            metadata[strainTrimmed]["5.rMLST"]["rMLSTSequenceType"] = reference
    return metadata


def dictionaryPreparer(genomes):
    """Populates dictionaries with 'name' and a 0 value to be incremented in the rMLST parsing function"""
    global genomeCovered, testDict
    for genome in genomes:
        genomeCovered[genome] = 0
        testDict[genome] = 0


def rMLSTsizer(metadata, sampleNames):
    """Adds the data stored in testDict (number of alleles with a sequence match to the profile),
    and genomeCovered (total number of identical nucleotides)"""
    for name in sampleNames:
        metadata[name]["1.General"]["CoreGeneAnalysisLength"] = "%s alleles (%s total nucleotides) matched at " \
                                                                "100%% identity" % (testDict[name], genomeCovered[name])
    return metadata


def functionsGoNOW(sampleNames, path, date, metadata, refFilesPath):
    """Commenting is subpar in this script, as I am using code written by Mike,
    so I can't really comment on it very well"""
    '''Yeah, well if I didn't have to change it so much to work it would have been commented better Adam'''
    print "\nPerforming rMLST analyses."
    rMLSTgenes = refFilesPath + "/rMLST/alleles/"
    make_path("%s/tmp" % path)
    print "\nFinding rMLST alleles."
    dictionaryPreparer(sampleNames)
    plusdict = blaster(rMLSTgenes, sampleNames, path, path, date, refFilesPath)
    additionalMetadata = determineReferenceGenome(plusdict, path, metadata, refFilesPath)
    moreMetadata = determineSubtype(plusdict, path, additionalMetadata, refFilesPath)
    allMetadata = rMLSTsizer(moreMetadata, sampleNames)
    # print json.dumps(genomeCovered, sort_keys=True, indent=4, separators=(',', ': '))
    # print json.dumps(testDict, sort_keys=True, indent=4, separators=(',', ': '))
    jsonReportR.jsonR(sampleNames, path, allMetadata, "Collection")
    return allMetadata
