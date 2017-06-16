#!/usr/bin/env python
from GeneSeekr import *

__author__ = 'adamkoziol'


class Univec(GeneSeekr):

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            genome = os.path.split(assembly)[1].split('.')[0]
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            make_path(sample[self.analysistype].reportdir)
            size = 0
            try:
                report = glob('{}{}*rawresults*'.format(sample[self.analysistype].reportdir, genome))[0]
                size = os.path.getsize(report)
            except IndexError:

                report = '{}{}_rawresults_{:}.csv'.format(sample[self.analysistype].reportdir, genome,
                                                          time.strftime("%Y.%m.%d.%H.%M.%S"))
            db = target.split('.')[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly,
                                           db=db,
                                           reward=1,
                                           penalty=-5,
                                           gapopen=3,
                                           gapextend=3,
                                           dust="yes",
                                           soft_masking="true",
                                           evalue=0.1,
                                           num_alignments=1000000,
                                           num_threads=24,
                                           outfmt="'6 qseqid sacc stitle positive mismatch gaps "
                                                  "evalue bitscore slen length'",
                                           out=report)
            # Save the blast command in the metadata
            sample[self.analysistype].blastcommand = str(blastn)
            if not os.path.isfile(report) or size == 0:
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
        blastdict = DictReader(open(report), fieldnames=self.univecnames, dialect='excel-tab')
        resultdict = {}
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            # noinspection PyTypeChecker
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            # Find the allele number and the text before the number for different formats
            target = '{},{},{}'.format(row['query_id'], row['subject_accession'], row['subject_title'])
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                resultdict.update({target: '{},{}'.format(str(percentidentity), row['positives'])})
            sample[self.analysistype].blastresults = resultdict

    def reporter(self):
        import operator
        combinedrow = ''
        for sample in self.metadata:
            row = ''
            # Populate the header with the appropriate data, including all the genes in the list of targets
            row += 'Contig,Accession,Title,PercentIdentity,Length\n'
            try:
                # Sort the results stored in the object and pull out the best match
                bestmatch = sorted(sample[self.analysistype].blastresults.items(),
                                   key=operator.itemgetter(1), reverse=True)
                # Iterate through all the matches in the dictionary
                for match in bestmatch:
                    # Update the string with the values in each match
                    row += '{},{}'.format(match[0], match[1])
                    row += '\n'
            except (AttributeError, KeyError):
                row = ''
            combinedrow += row
            # If the script is being run as part of the assembly pipeline, make a report for each sample
            if self.pipeline:
                if sample.general.bestassemblyfile != 'NA':
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
        self.univecnames = ['query_id', 'subject_accession', 'subject_title', 'positives', 'mismatches', 'gaps',
                            'evalue',  'bit_score', 'subject_length', 'alignment_length']
        GeneSeekr.__init__(self, inputobject)
