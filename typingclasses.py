#!/usr/bin/env python 3
from genesippr.genesippr import GeneSippr
from spadespipeline.GeneSeekr import *
__author__ = 'adamkoziol'


class Plasmids(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        data = 'Strain,Gene,PercentIdentity,Length,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    multiple = False
                    for name, identity in sample[self.analysistype].results.items():
                        if not multiple:
                            data += '{},{},{},{}\n'.format(name, identity,
                                                           len(sample[self.analysistype].sequences[name]),
                                                           sample[self.analysistype].avgdepth[name])
                        else:
                            data += ',{},{},{},{}\n'.format(name, identity,
                                                            len(sample[self.analysistype].sequences[name]),
                                                            sample[self.analysistype].avgdepth[name])
                        multiple = True
                else:
                    data += '\n'
            report.write(data)


class Resistance(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create a set of all the gene names without alleles or accessions e.g. sul1_18_AY260546 becomes sul1
        genedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(self.targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            for line in notes:
                # Ignore comment lines - they will break the parsing
                if line.startswith('#'):
                    continue
                # Split the line on colons e.g. QnrB53:Quinolone resistance: has three variables after the split:
                # gene(QnrB53), resistance(Quinolone resistance), and _(\n)
                gene, resistance, _ = line.split(':')
                # Set up the resistance dictionary
                genedict[gene] = resistance
        # Find unique gene names with the highest percent identity
        for sample in self.runmetadata.samples:
            if sample[self.analysistype].results:
                # Initialise a dictionary to store the unique genes, and their percent identities
                sample[self.analysistype].uniquegenes = dict()
                for name, identity in sample[self.analysistype].results.items():
                    # Split the name of the gene from the string e.g. ARR-2_1_HQ141279 yields ARR-2
                    genename = name.split('_')[0]
                    # Set the best observed percent identity for each unique gene
                    try:
                        # Pull the previous best identity from the dictionary
                        bestidentity = sample[self.analysistype].uniquegenes[genename]
                        # If the current identity is better than the old identity, save it
                        if float(identity) > float(bestidentity):
                            sample[self.analysistype].uniquegenes[genename] = float(identity)
                    # Initialise the dictionary if necessary
                    except KeyError:
                        sample[self.analysistype].uniquegenes[genename] = float(identity)
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        header = 'Strain,Resistance,Gene,Allele,Accession,PercentIdentity,Length,FoldCoverage\n'
        data = str()
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    # If there are multiple results for a sample, don't write the sample name in each line of the report
                    multiple = False
                    for name, identity in sample[self.analysistype].results.items():
                        # Split the name on underscores: ARR-2_1_HQ141279; gene: ARR-2, allele: 1, accession: HQ141279
                        genename, allele, accession = name.split('_')
                        # Retrieve the best identity for each gene
                        percentid = sample[self.analysistype].uniquegenes[genename]
                        # If the percent identity of the current gene matches the best percent identity, add it to
                        # the report - there can be multiple occurrences of genes e.g.
                        # sul1,1,AY224185,100.00,840 and sul1,2,CP002151,100.00,927 are both included because they
                        # have the same 100% percent identity
                        if float(identity) == percentid:
                            # Treat the initial vs subsequent results for each sample slightly differently - instead
                            # of including the sample name, use an empty cell instead
                            if not multiple:
                                # Populate the results
                                data += '{},{},{},{},{},{},{}\n'.format(
                                    genedict[genename],
                                    genename,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                            else:
                                data += ',{},{},{},{},{},{},{}\n'.format(
                                    genedict[genename],
                                    genename,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                            multiple = True
                else:
                    data += '\n'
            # Write the strings to the file
            report.write(header)
            report.write(data)


class Prophages(GeneSeekr):

    def reporter(self):
        with open(os.path.join(self.reportpath, 'prophages.csv'), 'w') as report:
            data = 'Strain,Gene,Host,PercentIdentity,PercentCovered,Contig,Location\n'
            # Set the required variables to load prophage data from a summary file
            targetpath = os.path.join(self.referencefilepath, self.analysistype)
            overview = glob(os.path.join(targetpath, '*.txt'))[0]
            fieldnames = ['id_prophage', 'file_name', 'host', 'host_id', 'number_of_prophages_in_host',
                          'start_position_of_prophage', 'end_position_of_prophage', 'length_of_prophage']
            for sample in self.metadata:
                # Create a set to ensure that genes are only entered into the report once
                genes = set()
                if sample.general.bestassemblyfile != 'NA':
                    # Open the prophage file as a dict - I do this here, as if I open it earlier, it looks like the
                    # file remains partially-read through for the next iteration. Something like prophagedata.seek(0)
                    # would probably work, but Dictreader objects don't have a .seek attribute
                    prophagedata = DictReader(open(overview), fieldnames=fieldnames, dialect='excel-tab')
                    if sample[self.analysistype].blastresults:
                        data += '{},'.format(sample.name)
                        # Allow for formatting multiple hits for the same sample
                        multiple = False
                        for result in sample[self.analysistype].blastresults:
                            gene = result['subject_id']
                            if gene not in genes:
                                if multiple:
                                    data += ','
                                # Iterate through the phage data in the dictionary
                                for phage in prophagedata:
                                    if phage['id_prophage'] == gene:
                                        # Add the data to the row
                                        data += '{},{},{},{},{},{}..{}\n' \
                                            .format(gene,
                                                    phage['host'],
                                                    result['percentidentity'],
                                                    result['alignment_fraction'] if float(result['alignment_fraction'])
                                                    <= 100 else '100.0',
                                                    result['query_id'],
                                                    result['low'],
                                                    result['high'])
                                genes.add(gene)
                                # Set multiple to true for any additional hits for this sample
                                multiple = True
                    else:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(data)


class Univec(GeneSeekr):

    def reporter(self):
        from Bio import SeqIO
        import re
        with open(os.path.join(self.reportpath, 'univec.csv'), 'w') as report:
            data = 'Strain,Gene,Description,PercentIdentity,PercentCovered,Contig,Location\n'
            for sample in self.metadata:
                if sample.general.bestassemblyfile != 'NA':
                    # Create a set to ensure that genes are only entered into the report once
                    genes = set()
                    if sample[self.analysistype].blastresults:
                        data += '{},'.format(sample.name)
                        # If multiple hits are returned for a sample, don't re-add the sample name on the next row
                        multiple = False
                        for result in sample[self.analysistype].blastresults:
                            gene = result['subject_id']
                            # Parse the reference file in order to extract the description of the BLAST hits
                            for entry in SeqIO.parse(sample[self.analysistype].combinedtargets, 'fasta'):
                                # Find the corresponding entry for the gene
                                if entry.id == gene:
                                    # Cut out the description from the entry.description using regex
                                    # e.g. for 'gnl|uv|X66730.1:1-2687-49 B.bronchiseptica plasmid pBBR1 genes for
                                    # mobilization and replication' only save the string after '2687-49'
                                    description = re.findall('\d+-\d+\s(.+)', entry.description)[0]
                                    # Don't add the same gene more than once to the report
                                    if gene not in genes:
                                        if multiple:
                                            data += ','
                                        data += '{},{},{},{},{},{}..{}\n' \
                                            .format(gene.split('|')[-1],
                                                    description,
                                                    result['percentidentity'],
                                                    result['alignment_fraction'] if float(result['alignment_fraction'])
                                                    <= 100 else '100.0',
                                                    result['query_id'],
                                                    result['low'],
                                                    result['high'])
                                        # Allow for the proper formatting
                                        multiple = True
                                        genes.add(gene)
                    else:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(data)


class Virulence(GeneSippr):

    def reporter(self):
        """
        Creates a report of the results
        """
        # Create a set of all the gene names without alleles or accessions e.g. sul1_18_AY260546 becomes sul1
        genedict = dict()
        # Load the notes file to a dictionary
        notefile = os.path.join(self.targetpath, 'notes.txt')
        with open(notefile, 'r') as notes:
            for line in notes:
                # Ignore comment lines - they will break the parsing
                if line.startswith('#'):
                    continue
                # Split the line on colons e.g. stx1Aa:  Shiga toxin 1, subunit A, variant a: has three variables after
                # the split: gene(stx1Aa), description(Shiga toxin 1, subunit A, variant a), and _(\n)
                try:
                    gene, description, _ = line.split(':')
                # There are exceptions to the parsing. Some lines only have one :, while others have three. Allow for
                # these possibilities.
                except ValueError:
                    try:
                        gene, description = line.split(':')
                    except ValueError:
                        gene, description, _, _ = line.split(':')
                # Set up the description dictionary
                genedict[gene] = description.replace(', ', '_').strip()
        # Find unique gene names with the highest percent identity
        for sample in self.runmetadata.samples:
            if sample[self.analysistype].results:
                # Initialise a dictionary to store the unique genes, and their percent identities
                sample[self.analysistype].uniquegenes = dict()
                for name, identity in sample[self.analysistype].results.items():
                    # Split the name of the gene from the string e.g. stx1:11:Z36899:11 yields stx1
                    genename = name.split(':')[0]
                    # Only allow matches of 100% identity for stx genes
                    if 'stx' in genename and float(identity) < 100.0:
                        pass
                    else:
                        # Set the best observed percent identity for each unique gene
                        try:
                            # Pull the previous best identity from the dictionary
                            bestidentity = sample[self.analysistype].uniquegenes[genename]
                            # If the current identity is better than the old identity, save it
                            if float(identity) > float(bestidentity):
                                sample[self.analysistype].uniquegenes[genename] = float(identity)
                        # Initialise the dictionary if necessary
                        except KeyError:
                            sample[self.analysistype].uniquegenes[genename] = float(identity)
                        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        data = 'Strain,Gene,Subtype/Allele,Description,Accession,PercentIdentity,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                if sample[self.analysistype].results:
                    # If there are multiple results for a sample, don't write the sample name in each line of the report
                    multiple = False
                    for name, identity in sorted(sample[self.analysistype].results.items()):
                        try:
                            # Split the name on colons: stx2A:63:AF500190:d; gene: stx2A, allele: 63, accession:
                            # AF500190, subtype: d
                            genename, allele, accession, subtype = name.split(':')
                        # Treat samples without a subtype e.g. icaC:intercellular adhesion protein C: differently.
                        # Extract the allele as the 'subtype', and the gene name, and accession as above
                        except ValueError:
                            genename, subtype, accession = name.split(':')
                        # Retrieve the best identity for each gene
                        percentid = sample[self.analysistype].uniquegenes[genename]
                        # If the percent identity of the current gene matches the best percent identity, add it to
                        # the report - there can be multiple occurrences of genes e.g.
                        # sul1,1,AY224185,100.00,840 and sul1,2,CP002151,100.00,927 are both included because they
                        # have the same 100% percent identity
                        if float(identity) == percentid:
                            # Treat the initial vs subsequent results for each sample slightly differently - instead
                            # of including the sample name, use an empty cell instead
                            if multiple:
                                data += ','
                            try:
                                description = genedict[genename]
                            except KeyError:
                                description = 'na'
                            # Populate the results
                            data += '{},{},{},{},{},{}\n'.format(
                                genename,
                                subtype,
                                description,
                                accession,
                                identity,
                                sample[self.analysistype].avgdepth[name])
                            multiple = True
                else:
                    data += '\n'
            # Write the strings to the file
            report.write(data)
