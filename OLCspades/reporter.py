#!/usr/bin/env python
from accessoryFunctions import *
__author__ = 'adamkoziol'


class Reporter(object):

    def reporter(self):
        from collections import OrderedDict
        printtime('Creating summary report', self.starttime)
        row = ''
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                data = OrderedDict([
                    ('SampleName', sample.name),
                    ('N50', sample.quast.N50.split('Count')[0]),
                    ('NumContigs', sample.mapping.Contigs),
                    ('TotalLength', sample.mapping.Bases.split('bp')[0]),
                    ('MeanInsertSize', sample.mapping.MeanInsertSize),
                    ('AverageCoverageDepth', sample.mapping.MeanCoveragedata.split("X")[0]),
                    ('ReferenceGenome', sample.rmlst.referencegenome),
                    ('RefGenomeAlleleMatches', str(sample.rmlst.matchestoreferencegenome)),
                    ('16sPhylogeny', sample['16s'].taxonomy),
                    ('rMLSTsequenceType', sample.rmlst.sequencetype),
                    ('MLSTsequencetype', sample.mlst.sequencetype),
                    ('MLSTmatches', str(sample.mlst.matchestosequencetype)),
                    ('coregenome', '{}/{}'.format(sample.coregenome.targetspresent, sample.coregenome.totaltargets)),
                    ('Serotype', sample.serotype.serotype),
                    ('geneSeekrProfile', ';'.join(result for result in sample.geneseekr.blastresults)
                        if sample.geneseekr.blastresults != 'NA' else 'NA'),
                    ('vtyperProfile', sample.vtyper.toxinprofile),
                    ('percentGC', sample.mapping.GcPercentage),
                    ('TotalPredictedGenes', str(sample.prodigal.predictedgenestotal)),
                    ('predictedgenesover3000bp', str(sample.prodigal.predictedgenesover3000bp)),
                    ('predictedgenesover1000bp', str(sample.prodigal.predictedgenesover1000bp)),
                    ('predictedgenesover500bp', str(sample.prodigal.predictedgenesover500bp)),
                    ('predictedgenesunder500bp', str(sample.prodigal.predictedgenesunder500bp)),
                    ('SequencingDate', sample.run.Date),
                    ('Investigator', sample.run.InvestigatorName),
                    ('TotalClustersinRun', str(sample.run.TotalClustersinRun)),
                    ('NumberofClustersPF', str(sample.run.NumberofClustersPF)),
                    ('PercentOfClusters', str(sample.run.PercentOfClusters)),
                    ('LengthofForwardRead', str(sample.run.forwardlength)),
                    ('LengthofReverseRead', str(sample.run.reverselength)),
                    ('Project', str(sample.run.SampleProject)),
                    ('PipelineVersion', self.commit)
                ])

                if not row:
                    row += ','.join([key for key, value in data.items()])
                row += '\n'
                row += ','.join(value for key, value in data.items())
        with open('{}/combinedMetadata.csv'.format(self.reportpath), 'wb') as metadatareport:
            metadatareport.write(row)

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.commit = inputobject.commit
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        self.reporter()
