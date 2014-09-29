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

count = 0


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

# Declare queues, list, and dict
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
# blastexist = {}

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
            # plusdict[path.group(3)] = {gene.group(1): 0}
            # print path.group(3), gene.group(1)
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
                                        # plusdict[genome][gene] = col
                                        plusdict[genome][gene[-2:]] = col
                                    # print genome, gene[-2:], col
                            threadlock.release()  # precaution for populate dictionary with GIL
            else:
                for genome in genomes:
                    for gene in genomes[genome]:
                        if genome not in plusdict:
                            plusdict[genome] = defaultdict(str)
                        if gene[-2:] not in plusdict[genome]:
                            # print gene, genome
                            # plusdict[genome][gene] = col
                            plusdict[genome][gene[-2:]] = "N"
                        # print genome, gene[-2:], col
                    # threadlock.release()  # precaution for populate dictionary with GIL
            dotter()
            mm.close()

            self.parsequeue.task_done()

def parsethreader(blastpath, genomes):
    global plusdict
    dotter()
    for i in range(len(genomes)):
        threads = blastparser(parsequeue)
        threads.setDaemon(True)
        threads.start()
    progress = len(blastpath)
    for xml in blastpath:
        # print xml
        handle = open(xml, 'r')
        mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        # time.sleep(0.05) # Previously used to combat open file error
        handle.close()
        parsequeue.put((xml, blastpath[xml], mm, progress))
        parsequeue.join()


def blaster(markers, strains, path, out, experimentName):
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
        genomes.append(genomeFile[0])
    # print genomes
    # if os.path.isdir(strains):
    #     genomes = glob.glob(strains + "*.fa")
    # elif os.path.isfile(strains):
    #     genomes = strains
    # else:
    #     print "The variable \"--genomes\" is not a folder or file"
    #     return
    sys.stdout.write("[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S")))
    # #push markers to threads
    makedbthreads(fastas)
    print "\n[%s] BLAST database(s) created" % (time.strftime("%H:%M:%S"))
    if os.path.isfile('%s/blastxmldict.json' % path):
        # print "[%s] Loading BLAST data from file" % (time.strftime("%H:%M:%S"))
        # sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
        blastpath = json.load(open('%s/blastxmldict.json' % path))
    else:
        # print "[%s] Now performing BLAST database searches" % (time.strftime("%H:%M:%S"))
        # sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
        # make blastn threads and retrieve xml file locations
        blastnthreads(fastas, genomes)
        json.dump(blastpath, open('%s/blastxmldict.json' % path, 'wb'), sort_keys=True, indent=4, separators=(',', ': '))
    print "\n[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    sys.stdout.write('[%s]' % (time.strftime("%H:%M:%S")))
    # for bpath in blastpath:
        # print fastas
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
            # for plusrow in plusdict[genomerow][generow]:
            row += ',' + str(plusdict[genomerow][generow])
    with open("%s/%s_results_%s.csv" % (out, experimentName, time.strftime("%H:%M:%S")), 'wb') as csvfile:
    # with open("%s/%s_results_%s.csv" % (out, experimentName, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        csvfile.write(csvheader)
        csvfile.write(row)
    print "\n[%s] Now parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    # print json.dumps(plusdict, sort_keys=True, indent=4)
    return plusdict


header = []
body = []

def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the dictionary responsible for storing the report data
sequenceTypes = defaultdict(make_dict)
strainTypes = defaultdict(make_dict)
mismatch = defaultdict(make_dict)


def determineSubtype(plusdict, path, metadata):
    # print json.dumps(plusdict, sort_keys=True, indent=4)
    with open("%s/referenceGenomes/referenceGenome_rMLSTprofiles.json" % path, "r") as referenceFile:
        sequenceTypes = json.load(referenceFile)
        referenceFile.close()

    # print json.dumps(referenceGenomeProfile, sort_keys=True, indent=4)
    # for strain in sorted(sequenceTypes):
    #     for gene, allele in sorted(sequenceTypes[strain].iteritems()):
    #         print strain, gene, allele

    # print json.dumps(plusdict, sort_keys=True, indent=4)
    #
    print "\nDetermining sequence types."
    for genome in plusdict:
        dotter()
        bestCount = 0
        for sType in sorted(sequenceTypes):
            count = 0
            for bactNum, allele in sorted(plusdict[genome].iteritems()):
            # for bactNum in sorted(plusdict[genome]):
                fullBact = "BACT0000%s" % bactNum
                if str(fullBact) in sequenceTypes[sType]:
                    # print genome, sType, fullBact, sequenceTypes[sType][str(fullBact)]
                    referenceAllele = sequenceTypes[sType][fullBact]
                    if re.search(" ", referenceAllele):
                        splitAllele = referenceAllele.split(" ")
                        if allele == splitAllele[0] or allele == splitAllele[1]:
                            count += 1
                        else:
                            mismatch[fullBact][allele] = sequenceTypes[sType][fullBact]
                    else:
                        if allele == referenceAllele:
                            count += 1
                        else:
                            mismatch[fullBact][allele] = sequenceTypes[sType][fullBact]
                else:
                    # print genome, sType, fullBact, "N"
                    mismatch[fullBact][allele] = "N"
                # if allele == sequenceTypes[sType][fullBact]:
                    # count += 1
                # else:
                    # mismatch[fullBact][allele] = sequenceTypes[sType][fullBact]
            if count > bestCount:
                bestCount = count
                strainTypes[genome].clear()
                # print genome, sType, json.dumps(mismatch, sort_keys=True, indent=4)
                if mismatch:
                    for geneName in sorted(mismatch):
                        for observedAllele, refAllele in sorted(mismatch[geneName].iteritems()):
                            strainTypes[genome][sType][bestCount][geneName][observedAllele] = refAllele
                else:
                    strainTypes[genome][sType] = bestCount
                    # print genome, sType, "no mismatches"
                # print genome, sType, count
            mismatch.clear()
    # print json.dumps(strainTypes, sort_keys=True, indent=4)
    for strain in sorted(strainTypes):
        strainTrimmed = re.split("_filteredAssembled", strain)[0]
        make_path("%s/%s/referenceGenome" % (path, strainTrimmed))
        for reference in sorted(strainTypes[strain]):
            shutil.copy("%s/referenceGenomes/%s.fasta" % (path, reference), "%s/%s/referenceGenome" % (path, strainTrimmed))
            metadata[strainTrimmed]["1.General"]["rMLSTmatchestoRef"] = strainTypes[strain][reference]
            # print strain, reference


def functionsGoNOW(sampleNames, path, date, metadata):
    print "\nPerforming rMLST analyses."
    rMLSTgenes = path + "/rMLST/alleles/"
    plusdict = blaster(rMLSTgenes, sampleNames, path, path, date)
    determineSubtype(plusdict, path, metadata)
    # print rMLSTgenes
    # print path, sampleNames