#!/usr/bin/env python
from GeneSeekr import *
__author__ = 'adamkoziol'


class ResFinder(GeneSeekr):

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            make_path(sample[self.analysistype].reportdir)
            try:
                report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
                size = os.path.getsize(report)
                if size == 0:
                    os.remove(report)
                    report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                              time.strftime("%Y.%m.%d.%H.%M.%S"))
            except IndexError:
                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
            db = target.split('.')[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly, db=db, evalue='1E-20', num_alignments=1000000,
                                           num_threads=12,
                                           outfmt='"6 qseqid sseqid positive mismatch gaps '
                                                  'evalue bitscore slen length qstart qend"',
                                           out=report)
            if not os.path.isfile(report):
                try:
                    blastn()
                except:
                    self.blastqueue.task_done()
                    self.blastqueue.join()
                    try:
                        os.remove(report)
                    except IOError:
                        pass
                    raise
            # Run the blast parsing module
            self.blastparser(report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.resfinderfields, dialect='excel-tab')
        resultdict = {}
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            target = row['subject_id']
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Update the dictionary with the target, query, and percent identity
                resultdict.update({target: {row['query_id']: percentidentity}})
            sample[self.analysistype].blastresults = resultdict
        if not resultdict:
            sample[self.analysistype].blastresults = 'NA'

    def csvwriter(self):
        combinedrow = ''
        for sample in self.metadata:
            row = ''
            if sample[self.analysistype].targetnames != 'NA':
                try:
                    resistance = ''
                    notes = open('{}notes.txt'.format(sample[self.analysistype].targetpath), 'rb').readlines()
                    # Populate the header with the appropriate data, including all the genes in the list of targets
                    row += 'Contig,Accession,Gene,Resistance,PercentIdentity\n'
                    # Create a set of all the genes present in the blast results - list comprehension, and
                    # split on underscore
                    geneset = set([accession.split("_")[0] for accession in sample[self.analysistype].blastresults])
                    # Iterate through the list of genes
                    for gene in sorted(geneset):
                        # Find the resistance associated with the gene in the notes file
                        for line in notes:
                            # If the gene matches the gene in the line of interest
                            if gene == line.split(':')[0]:
                                # Set :resistance appropriately
                                resistance = line.split(':')[-2].replace(' resistance', '')
                        # Determine the best identity for all the results for the gene: max() of a list comprehension
                        # that returns all the %ids in the blast results for that particular gene
                        maxpercentid = max([percentid.items()[0][1] for accession, percentid in
                                            sample[self.analysistype].blastresults.items() if gene in accession])
                        # Iterate through all the blast results
                        for accession in sample[self.analysistype].blastresults:
                            # If the %id and the gene name match the gene of interest and the max %id
                            if sample[self.analysistype].blastresults[accession].items()[0][1] == maxpercentid and \
                                            accession.split("_")[0] == gene:
                                # Add the appropriate values to the row based on the header
                                row += ','.join([sample[self.analysistype].blastresults[accession].items()[0][0],
                                                 accession.split("_")[-1], gene, resistance, str(maxpercentid)])
                                row += '\n'
                # If there are no results, row stays empty
                except AttributeError:
                    row = ''
            combinedrow += row
            # If the script is being run as part of the assembly pipeline, make a report for each sample
            if self.pipeline:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[self.analysistype].reportdir != 'NA':
                        # Open the report
                        with open('{}{}_{}.csv'.format(sample[self.analysistype].reportdir, sample.name,
                                                       self.analysistype), 'wb') as report:
                            # Write the row to the report
                            report.write(row)
            try:
                # Remove the messy blast results from the object
                delattr(sample[self.analysistype], "blastresults")
            except KeyError:
                pass
        # Create the report containing all the data from all samples
        with open('{}{}.csv'.format(self.reportpath, self.analysistype), 'wb') \
                as combinedreport:
            combinedreport.write(combinedrow)

    def __init__(self, inputobject):
        # qseqid sacc stitle positive mismatch gaps evalue bitscore slen length
        self.resfinderfields = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps', 'evalue', 'bit_score',
                                'subject_length', 'alignment_length', 'query_start', 'query_end']
        GeneSeekr.__init__(self, inputobject)
