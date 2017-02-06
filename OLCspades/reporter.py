#!/usr/bin/env python
from accessoryFunctions import *
__author__ = 'adamkoziol'


class Reporter(object):

    def reporter(self):
        from collections import OrderedDict
        printtime('Creating summary report', self.starttime)
        row = ''
        # Create a dictionary of tuples to be printed in the final report
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
                    ('16sPhylogeny', sample['sixteenS'].taxonomy),
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
            else:
                if not row:
                    data = ['SampleName', 'N50', 'NumContigs', 'TotalLength', 'MeanInsertSize', 'AverageCoverageDepth',
                            'ReferenceGenome', 'RefGenomeAlleleMatches', '16sPhylogeny', 'rMLSTsequenceType',
                            'MLSTsequencetype', 'MLSTmatches', 'coregenome', 'Serotype', 'geneSeekrProfile',
                            'vtyperProfile', 'percentGC', 'TotalPredictedGenes', 'predictedgenesover3000bp',
                            'predictedgenesover1000bp', 'predictedgenesover500bp', 'predictedgenesunder500bp',
                            'SequencingDate', 'Investigator', 'TotalClustersinRun', 'NumberofClustersPF',
                            'PercentOfClusters', 'LengthofForwardRead', 'LengthofReverseRead', 'Project',
                            'PipelineVersion']
                    row += ','.join(data)
                row += '\n{}'.format(sample.name)
        with open('{}/combinedMetadata.csv'.format(self.reportpath), 'wb') as metadatareport:
            metadatareport.write(row)

    def database(self):
        """
        Enters all the metadata into a database
        """
        import sqlite3
        try:
            os.remove('{}/metadatabase.sqlite'.format(self.reportpath))
        except OSError:
            pass
        # Set the name of the database
        db = sqlite3.connect('{}/metadatabase.sqlite'.format(self.reportpath))
        # Create a cursor to allow access to the database
        cursor = db.cursor()
        # Set up the db
        cursor.execute('''
          CREATE TABLE IF NOT EXISTS Samples (
            id     INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
            name   TEXT UNIQUE
          )

        ''')
        # Create a variable to store the names of the header values for each individual table
        # This will store a set of all the headers from all the strains, as there can be some variability present, as
        # not all analyses are available for all taxonomic groups
        columns = dict()
        for sample in self.metadata:
            # Insert each strain name into the Samples table
            cursor.execute('''
              INSERT OR IGNORE INTO Samples (name)
              VALUES ( ? )
            ''', (sample.name, ))
            # Each header in the .json file represents a major category e.g. ARMI, geneseekr, commands, etc. and
            # will be made into a separate table
            for header in sample.datastore.items():
                # Set the table name
                tablename = header[0].replace('.', '_')
                # Allow for certain analyses, such as core genome, not being performed on all strains
                try:
                    # Create the table (if it doesn't already exist)
                    cursor.execute('''
                      CREATE TABLE IF NOT EXISTS {} (
                        sample_id INTEGER
                      )
                      '''.format(tablename))
                    # Key and value: data description and data value e.g. targets present: 1012, etc.
                    for key, value in header[1].datastore.items():
                        # Add the data header to the dictionary
                        # Clean the column names so there are no issues entering names into the database
                        cleanedcolumn = self.columnclean(key)
                        try:
                            columns[tablename].add(str(cleanedcolumn))
                        # Initialise the dictionary the first time a table name is encountered
                        except KeyError:
                            columns[tablename] = set()
                            columns[tablename].add(str(cleanedcolumn))
                except (AttributeError, IndexError):
                    pass
        # Iterate through the dictionary containing all the data headers
        for table, setofheaders in columns.items():
            # Each header will be used as a column in the appropriate table
            for cleanedcolumn in setofheaders:
                # Alter the table by adding each header as a column
                cursor.execute('''
                  ALTER TABLE {}
                  ADD COLUMN {} TEXT
                '''.format(table, cleanedcolumn))
            # Iterate through the samples and pull out the data for each table/column
            for sample in self.metadata:
                # Find the id associated with each sample in the Sample table
                cursor.execute('''
                  SELECT id from Samples WHERE name=?
                ''', (sample.name,))
                sampleid = cursor.fetchone()[0]
                # Add the sample_id to the table
                cursor.execute('''
                  INSERT OR IGNORE INTO {}
                  (sample_id) VALUES ("{}")
                 '''.format(table, sampleid))
                # Add the data to the table
                try:
                    # Find the data for each table/column
                    for item in sample[table].datastore.items():
                        # Clean the names
                        cleanedcolumn = self.columnclean(str(item[0]))
                        # Add the data to the column of the appropriate table,
                        # where the sample_id matches the current strain
                        cursor.execute('''
                          UPDATE {}
                          SET {} = ?
                          WHERE sample_id = {}
                          '''.format(table, cleanedcolumn, sampleid), (str(item[1]), ))
                except KeyError:
                    pass
        # Commit the changes to the database
        db.commit()

    @staticmethod
    def columnclean(column):
        """
        Modifies column header format to be importable into a database
        :param column: raw column header
        :return: cleanedcolumn: reformatted column header
        """
        cleanedcolumn = str(column) \
            .replace('%', 'percent') \
            .replace('(', '_') \
            .replace(')', '') \
            .replace('As', 'Adenosines') \
            .replace('Cs', 'Cytosines') \
            .replace('Gs', 'Guanines') \
            .replace('Ts', 'Thymines') \
            .replace('Ns', 'Unknowns') \
            .replace('index', 'adapterIndex')
        return cleanedcolumn

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.commit = inputobject.commit
        self.reportpath = inputobject.reportpath
        self.starttime = inputobject.starttime
        self.reporter()
        # Create a database to store all the metadata
        self.database()
