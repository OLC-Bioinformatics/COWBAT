#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import filer, GenObject, MetadataObject, printtime, make_path, write_to_logfile
from spadespipeline.GeneSeekr import GeneSeekr
from genesippr.genesippr import GeneSippr
import confindr.confindr as confinder
from biotools import bbtools
from csv import DictReader
from glob import glob
import pandas
import shutil
import os

__author__ = 'adamkoziol'


class Quality(object):

    def contamination_finder(self):
        """

        """
        printtime('Calculating contamination in reads', self.start)
        report = os.path.join(self.path, 'confinder', 'confindr_report.csv')
        if not os.path.isfile(report):
            # Create an object to store attributes to pass to confinder
            args = MetadataObject
            args.input_directory = self.path
            args.output_name = os.path.join(self.path, 'confinder')
            args.databases = os.path.join(self.reffilepath, 'ConFindr', 'databases')
            args.forward_id = '_R1'
            args.reverse_id = '_R2'
            args.threads = self.cpus
            args.kmer_size = 31
            args.number_subsamples = 5
            args.subsample_depth = 20
            args.kmer_cutoff = 2
            try:
                shutil.rmtree(args.output_name)
            except IOError:
                pass
            # Create a detector object
            confinder.run_mashsippr(args.input_directory,
                                    args.output_name,
                                    args.databases)
            # Open the output report file.
            with open(os.path.join(report), 'w') as f:
                f.write('Sample,Genus,NumContamSNVs,NumUniqueKmers,CrossContamination,ContamStatus\n')
            paired_reads = confinder.find_paired_reads(args.input_directory,
                                                       forward_id=args.forward_id,
                                                       reverse_id=args.reverse_id)
            # Perform contamination detection on each set of paired reads
            for pair in paired_reads:
                sample_name = os.path.basename(list(filer(pair))[0])
                printtime('Beginning analysis of sample {}...\n'.format(sample_name), self.start, '\033[1;34m')
                genus = confinder.read_mashsippr_output(os.path.join(args.output_name, 'reports', 'mash.csv'),
                                                        sample_name)
                confinder.find_contamination(pair, args, genus)
            printtime('Contamination detection complete!', self.start)
        # Load the confindr report into a dictionary using pandas
        # https://stackoverflow.com/questions/33620982/reading-csv-file-as-dictionary-using-pandas
        confindr_results = pandas.read_csv(report, index_col=0).T.to_dict()
        # Find the results for each of the samples
        for sample in self.metadata:
            # Create a GenObject to store the results
            sample.confinder = GenObject()
            # Iterate through the dictionary to find the outputs for each sample
            for line in confindr_results:
                # If the current line corresponds to the sample of interest
                if sample.name in line:
                    # Set the values using the appropriate keys as the attributes
                    sample.confinder.genus = confindr_results[line]['Genus']
                    sample.confinder.num_contaminated_snvs = confindr_results[line]['NumContamSNVs']
                    sample.confinder.unique_kmers = confindr_results[line]['NumUniqueKmers']
                    sample.confinder.cross_contamination = confindr_results[line]['CrossContamination']
                    sample.confinder.contam_status = confindr_results[line]['ContamStatus']
                    if sample.confinder.contam_status is True:
                        sample.confinder.contam_status = 'Contaminated'
                    elif sample.confinder.contam_status is False:
                        sample.confinder.contam_status = 'Clean'

    def estimate_genome_size(self):
        """
        Use kmercountexact from the bbmap suite of tools to estimate the size of the genome
        """
        printtime('Estimating genome size using kmercountexact', self.start)
        for sample in self.metadata:
            # Initialise the name of the output file
            sample[self.analysistype].peaksfile = os.path.join(sample[self.analysistype].outputdir, 'peaks.txt')
            # Run the kmer counting command
            out, err, cmd = bbtools.kmercountexact(forward_in=sorted(sample.general.fastqfiles)[0],
                                                   peaks=sample[self.analysistype].peaksfile,
                                                   returncmd=True,
                                                   threads=self.cpus)
            # Set the command in the object
            sample[self.analysistype].kmercountexactcmd = cmd
            # Extract the genome size from the peaks file
            sample[self.analysistype].genomesize = bbtools.genome_size(sample[self.analysistype].peaksfile)
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)

    def error_correction(self):
        """
        Use tadpole from the bbmap suite of tools to perform error correction of the reads
        """
        printtime('Error correcting reads', self.start)
        for sample in self.metadata:
            sample.general.trimmedcorrectedfastqfiles = [fastq.split('.fastq.gz')[0] + '_trimmed_corrected.fastq.gz'
                                                         for fastq in sorted(sample.general.fastqfiles)]
            out, err, cmd = bbtools.tadpole(forward_in=sorted(sample.general.trimmedfastqfiles)[0],
                                            forward_out=sample.general.trimmedcorrectedfastqfiles[0],
                                            returncmd=True,
                                            mode='correct',
                                            threads=self.cpus)
            # Set the command in the object
            sample[self.analysistype].errorcorrectcmd = cmd
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)

    def normalise_reads(self):
        """
        Use bbnorm from the bbmap suite of tools to perform read normalisation
        """
        printtime('Normalising reads to a kmer depth of 100', self.start)
        for sample in self.metadata:
            # Set the name of the normalised read files
            sample.general.normalisedreads = [fastq.split('.fastq.gz')[0] + '_normalised.fastq.gz'
                                              for fastq in sorted(sample.general.fastqfiles)]
            # Run the normalisation command
            out, err, cmd = bbtools.bbnorm(forward_in=sorted(sample.general.trimmedcorrectedfastqfiles)[0],
                                           forward_out=sample.general.normalisedreads[0],
                                           returncmd=True,
                                           threads=self.cpus)
            sample[self.analysistype].normalisecmd = cmd
            write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)

    def merge_pairs(self):
        """
        Use bbmerge from the bbmap suite of tools to merge paired-end reads
        """
        printtime('Merging paired reads', self.start)
        for sample in self.metadata:
            # Can only merge paired-end
            if len(sample.general.fastqfiles) == 2:
                # Set the name of the merged, and unmerged files
                sample.general.mergedreads = \
                    os.path.join(sample.general.outputdirectory, '{}_paired.fastq.gz'.format(sample.name))
                sample.general.unmergedforward = \
                    os.path.join(sample.general.outputdirectory, '{}_unpaired_R1.fastq.gz'.format(sample.name))
                sample.general.unmergedreverse = \
                    os.path.join(sample.general.outputdirectory, '{}_unpaired_R2.fastq.gz'.format(sample.name))
                # Run the merging command
                out, err, cmd = bbtools.bbmerge(forward_in=sample.general.normalisedreads[0],
                                                merged_reads=sample.general.mergedreads,
                                                returncmd=True,
                                                outu1=sample.general.unmergedforward,
                                                outu2=sample.general.unmergedreverse,
                                                threads=self.cpus)
                sample[self.analysistype].bbmergecmd = cmd
                write_to_logfile(out, err, self.logfile, sample.general.logout, sample.general.logerr, None, None)
            else:
                sample.general.mergedreads = sorted(sample.general.trimmedcorrectedfastqfiles)[0]

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.path = inputobject.path
        self.analysistype = 'quality'
        self.reffilepath = inputobject.reffilepath
        self.cpus = inputobject.cpus
        self.logfile = inputobject.logfile
        # Initialise the quality attribute in the metadata object
        for sample in self.metadata:
            setattr(sample, self.analysistype, GenObject())


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
                try:
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
                except KeyError:
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
            try:
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
            except KeyError:
                pass
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        header = 'Strain,Resistance,Gene,Allele,Accession,PercentIdentity,Length,FoldCoverage\n'
        data = str()
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                try:
                    if sample[self.analysistype].results:
                        # Create an attribute to store the string for the eventual pipeline report
                        sample[self.analysistype].pipelineresults = list()
                        # If there are multiple results for a sample, don't write the name in each line of the report
                        multiple = False
                        for name, identity in sample[self.analysistype].results.items():
                            # Split the name on '_'s: ARR-2_1_HQ141279; gene: ARR-2, allele: 1, accession: HQ141279
                            try:
                                genename, allele, accession = name.split('_')
                            # Some names have a slightly different naming scheme:
                            # tet(44)_1_NZ_ABDU01000081; gene: tet(44), allele: 1, accession: NZ_ABDU01000081
                            except ValueError:
                                genename, allele, preaccession, postaccession = name.split('_')
                                accession = '{preaccess}_{postaccess}'.format(preaccess=preaccession,
                                                                              postaccess=postaccession)
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
                                # Populate the results
                                data += '{},{},{},{},{},{},{}\n'.format(
                                    genedict[genename],
                                    genename,
                                    allele,
                                    accession,
                                    identity,
                                    len(sample[self.analysistype].sequences[name]),
                                    sample[self.analysistype].avgdepth[name])
                                sample[self.analysistype].pipelineresults.append(
                                    '{rgene} ({pid}%) {rclass}'.format(rgene=genename,
                                                                       pid=identity,
                                                                       rclass=genedict[genename])
                                )
                                multiple = True
                    else:
                        data += '\n'
                except KeyError:
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
                    try:
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
                                                        result['alignment_fraction'] if float(
                                                            result['alignment_fraction']) <= 100 else '100.0',
                                                        result['query_id'],
                                                        result['low'],
                                                        result['high'])
                                    genes.add(gene)
                                    # Set multiple to true for any additional hits for this sample
                                    multiple = True
                        else:
                            data += '{}\n'.format(sample.name)
                    except KeyError:
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
                    try:
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
                                                        result['alignment_fraction'] if float(
                                                            result['alignment_fraction'])
                                                        <= 100 else '100.0',
                                                        result['query_id'],
                                                        result['low'],
                                                        result['high'])
                                            # Allow for the proper formatting
                                            multiple = True
                                            genes.add(gene)
                        else:
                            data += '{}\n'.format(sample.name)
                    except KeyError:
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
            try:
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
            except KeyError:
                pass
        # Create the path in which the reports are stored
        make_path(self.reportpath)
        # Initialise strings to store the results
        data = 'Strain,Gene,Subtype/Allele,Description,Accession,PercentIdentity,FoldCoverage\n'
        with open(os.path.join(self.reportpath, self.analysistype + '.csv'), 'w') as report:
            for sample in self.runmetadata.samples:
                data += sample.name + ','
                try:
                    if sample[self.analysistype].results:
                        # If there are many results for a sample, don't write the sample name in each line of the report
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
                except KeyError:
                    data += '\n'
            # Write the strings to the file
            report.write(data)
