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
        os.chmod(inPath, 0777)
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
            global blastpath, plusdict
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


class blastparser(threading.Thread): # records, genomes):
    def __init__(self, parsequeue):
        self.parsequeue = parsequeue
        threading.Thread.__init__(self)
    def run(self):
        while True:
            global plusdict, genedict
            xml, genomes, mm, num = self.parsequeue.get()
            records = NCBIXML.parse(mm)
            numhsp = sum(line.count('<Hsp>') for line in iter(mm.readline, ""))
            if numhsp >= 1:
                mm.seek(0)
                for record in records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            threadlock.acquire()  # precaution
                            col = 'N'
                            if hsp.identities == alignment.length:
                                col = alignment.title.split('_')[-1]  # MLST type
                                col = re.sub(" No definition line", "", col)
                            for genome in genomes:
                                for gene in genomes[genome]:
                                    if genome not in plusdict:
                                        plusdict[genome] = defaultdict(str)
                                    if gene[-2:] not in plusdict[genome]:
                                        plusdict[genome][gene[-2:]] = col
                            threadlock.release()  # precaution for populate dictionary with GIL
            else:
                for genome in genomes:
                    for gene in genomes[genome]:
                        if genome not in plusdict:
                            plusdict[genome] = defaultdict(str)
                        if gene[-2:] not in plusdict[genome]:
                            plusdict[genome][gene[-2:]] = "N"
            dotter()
            mm.close()
            self.parsequeue.task_done()


def parsethreader(bPath, Genomes):
    """Sets up multithreaded blast"""
    global plusdict
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
    with open("%s/reports/%s_rMLSTresults.csv" % (out, experimentName), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    return plusdict


def determineReferenceGenome(plusdict, path, metadata, refFilesPath):
    """Find the reference genome with the closest number of identical alleles to each
    sequenced strain"""
    # Retrieve the rMLST profiles of the reference genomes
    with open("%s/referenceGenomes/referenceGenome_rMLSTprofiles.json" % refFilesPath, "r") as referenceFile:
        sequenceTypes = json.load(referenceFile)
        referenceFile.close()
    print "\nDetermining closest reference genome."
    for genome in plusdict:
        dotter()
        bestCount = 0
        for sType in sorted(sequenceTypes):
            count = 0
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
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
            shutil.copy("%s/referenceGenomes/%s.fasta" % (refFilesPath, reference), "%s/%s/referenceGenome" % (path, strainTrimmed))
            # print json.dumps(strainTypes[strain][reference], sort_keys=True, indent=4, separators=(',', ': '))
            metadata[strainTrimmed]["6.rMLSTmatchestoRef"] = strainTypes[strain][reference]
            # metadata[strainTrimmed]["6.rMLSTmatchestoRef"]["NumIdenticalAlleles"] = bestCount


    return metadata


def determineSubtype(plusdict, path, metadata, refFilesPath):
    """Same as above (determineReferenceGenome), but this determines the rMLST sequence type
    Comments are as above, as this is very similar"""
    profile = open("%s/rMLST/profile/rMLST_scheme.txt" % refFilesPath)
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
    # Populate the dictionary of sequenceType => geneName => alleleValue
    for line in body:
        for i in range(len(header) + 1):
            sequenceTypesMLST[int(line[0])][header[i - 1]] = line[i]
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
                else:
                    mismatchMLST[fullBact][allele] = sequenceTypesMLST[sType][fullBact]
            if count > bestCount:
                bestCount = count
                strainTypesMLST[genome].clear()
                if mismatchMLST:
                    for geneName in sorted(mismatchMLST):
                        for observedAllele, refAllele in sorted(mismatchMLST[geneName].iteritems()):
                            # strainTypesMLST[genome][sType]["NumIdenticalAlleles"][bestCount]["gene"][geneName]["observedAllele"][observedAllele]["referenceAllele"] = refAllele
                            strainTypesMLST[genome][str(sType)]["rMLSTIdenticalAlleles"] = bestCount
                            strainTypesMLST[genome][str(sType)]["rMLSTMismatchDetails"][geneName]["observedAllele"][observedAllele]["referenceAllele"] = refAllele
                            # print json.dumps(strainTypesMLST[genome][sType], sort_keys=True, indent=4, separators=(',', ': '))
                else:
                    strainTypesMLST[genome][str(sType)]["rMLSTIdenticalAlleles"] = bestCount
            mismatchMLST.clear()
    # Somehow, there are some files which do not have a single allele in common with any
    # of the profiles present in the rMLST scheme. These are addressed here.
    for genome in plusdict:
        # If the genome key is not in the strainTypesMLST dictionary
        if genome not in sorted(strainTypesMLST):
            strainTrimmed = re.split("_filteredAssembled", genome)[0]
            metadata[strainTrimmed]["5.rMLST"]["rMLSTSequenceType"] = "N/A"
    # Populates the metadata dictionary with data from strains present in the strainTypesMLST dictionary
    for strain in sorted(strainTypesMLST):
        strainTrimmed = re.split("_filteredAssembled", strain)[0]
        for INTreference in sorted(strainTypesMLST[strain]):
            reference = str(INTreference)
            # metadata[strainTrimmed]["1.General"]["rMLSTSequenceType"]["sequenceType"][reference] = strainTypesMLST[strain][reference]
            # print reference, strainTrimmed
            metadata[strainTrimmed]["5.rMLST"] = strainTypesMLST[strain][reference]
            # print json.dumps(metadata[strainTrimmed]["5.rMLST"], sort_keys=True, indent=4, separators=(',', ': '))
            # print json.dumps(strainTypesMLST[strain][reference], sort_keys=True, indent=4, separators=(',', ': '))
            metadata[strainTrimmed]["5.rMLST"]["rMLSTSequenceType"] = reference
    return metadata


def functionsGoNOW(sampleNames, path, date, metadata, refFilesPath):
    """Commenting is subpar in this script, as I am using code written by Mike,
    so I can't really comment on it very well"""
    '''Yeah, well if I didn't have to change it so much to work it would have been commented better Adam'''
    print "\nPerforming rMLST analyses."
    rMLSTgenes = refFilesPath + "/rMLST/alleles/"
    make_path("%s/tmp" % path)
    print "\nFinding rMLST alleles."
    plusdict = blaster(rMLSTgenes, sampleNames, path, path, date, refFilesPath)
    additionalMetadata = determineReferenceGenome(plusdict, path, metadata, refFilesPath)
    moreMetadata = determineSubtype(plusdict, path, additionalMetadata, refFilesPath)
    # print json.dumps(moreMetadata, sort_keys=True, indent=4, separators=(',', ': '))
    return moreMetadata
