#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import MetadataObject, GenObject, make_path
from spadespipeline import metadataReader
from argparse import ArgumentParser
import multiprocessing
from time import time
import pytest
import shutil
import sys
import os

testpath = os.path.abspath(os.path.dirname(__file__))
scriptpath = os.path.join(testpath, '..')
sys.path.append(scriptpath)
from assembly_pipeline import RunSpades

__author__ = 'adamkoziol'


@pytest.fixture()
def variables():
    v = ArgumentParser()
    v.path = os.path.join(testpath, 'testdata')
    v.sequencepath = v.path
    v.referencefilepath = os.path.join(v.path, 'databases')
    v.customsamplesheet = os.path.join(v.path, 'SampleSheet.csv')
    v.debug = True
    v.numreads = 2
    v.kmerrange = '21'
    v.preprocess = False
    v.basicassembly = True
    v.threads = multiprocessing.cpu_count()
    return v


@pytest.fixture()
def method_init(variables):
    method = RunSpades(variables, b'', time(), scriptpath)
    return method


@pytest.fixture()
def read_metadata(variables):
    metadata = metadataReader.MetadataReader(variables)
    return metadata


method = method_init(variables())


def test_sistr(variables):
    metadata = MetadataObject()
    method.runmetadata.samples = list()
    fasta = os.path.join(variables.path, 'NC_003198.fasta')
    metadata.name = os.path.split(fasta)[1].split('.')[0]
    # Initialise the general and run categories
    metadata.general = GenObject()
    metadata.run = GenObject()
    metadata.general.fastqfiles = list()
    # Set the destination folder
    outputdir = os.path.join(variables.sequencepath, metadata.name)
    make_path(outputdir)
    # Add the output directory to the metadata
    metadata.general.outputdirectory = outputdir
    metadata.run.outputdirectory = outputdir
    metadata.general.bestassemblyfile = True
    # Initialise an attribute to store commands
    metadata.commands = GenObject()
    # Assume that all samples are Salmonella
    metadata.general.referencegenus = 'Salmonella'
    # Set the .fasta file as the best assembly
    metadata.general.bestassemblyfile = fasta
    method.runmetadata.samples.append(metadata)
    method.sistr()
    for sample in method.runmetadata.samples:
        assert sample.sistr.cgmlst_genome_match == 'SAL_BA2732AA'
    variable_update()


def variable_update():
    global method
    method = method_init(variables())


def test_basic_link(variables):
    method.helper()
    assert os.path.islink(os.path.join(variables.path, 'NC_002695', 'NC_002695_R1.fastq.gz'))


def test_metadata():
    method.helper()
    for sample in method.runmetadata.samples:
        assert sample.name == 'NC_002695'


def test_basic_read_length():
    method.helper()
    for sample in method.runmetadata.samples:
        assert sample.run.forwardlength == 301


def test_quality_object():
    method.create_quality_object()


def test_raw_fastqc_paired():
    method.fastqc_raw()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'Raw', 'NC_002695_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_raw_fastqc_forward():
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'Raw', 'NC_002695_R1_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_quality_trim(variables):
    method.quality_trim()
    outfile = os.path.join(variables.path, 'NC_002695', 'NC_002695_R1_trimmed.fastq.gz')
    size = os.stat(outfile)
    assert size.st_size > 0


def test_trimmed_fastqc():
    method.fastqc_trimmed()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'Trimmed', 'NC_002695_R1_trimmed_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_error_correction(variables):
    method.error_correct()
    assert os.path.isfile(os.path.join(variables.path, 'NC_002695', 'NC_002695_R1_trimmed_corrected.fastq.gz'))


def test_confindr():
    method.contamination_detection()
    for sample in method.runmetadata.samples:
        assert sample.confinder.unique_kmers >= 575


def test_trimmed_corrected_fastqc():
    method.fastqc_trimmedcorrected()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'trimmedcorrected',
                               'NC_002695_R1_trimmed_corrected_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_read_normalisation(variables):
    method.normalise_reads()
    assert os.path.isfile(os.path.join(variables.path, 'NC_002695', 'NC_002695_R1_normalised.fastq.gz'))


def test_normalised_fastqc():
    method.fastqc_normalised()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'normalised',
                               'NC_002695_R1_normalised_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_read_merging(variables):
    method.merge_reads()
    assert os.path.isfile(os.path.join(variables.path, 'NC_002695', 'NC_002695_paired.fastq.gz'))


def test_merged_fastqc():
    method.fastqc_merged()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'fastqc', 'merged',
                               'NC_002695_paired_fastqc.zip')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_spades():
    method.run_spades()
    for sample in method.runmetadata.samples:
        outfile = os.path.join(sample.general.outputdirectory, 'spades_output', 'contigs.fasta')
        size = os.stat(outfile)
        assert size.st_size > 0


def test_qualimap():
    method.qualimap()
    for sample in method.runmetadata.samples:
        assert int(sample.mapping.Contigs) >= 500


def test_quality_features():
    method.quality_features()
    for sample in method.runmetadata.samples:
        assert sample.quality_features.n50 >= 845


def test_prodigal():
    method.prodigal()
    for sample in method.runmetadata.samples:
        assert sample.prodigal.predictedgenesover1000bp >= 65


def test_mash():
    method.mash()
    for sample in method.runmetadata.samples:
        assert sample.mash.closestrefseq == 'Escherichia coli O157:H7 str. Sakai'


def test_rmlst():
    method.rmlst()
    for sample in method.runmetadata.samples:
        assert sample.rmlst.sequencetype == '2124'


def test_sixteens():
    method.sixteens()
    for sample in method.runmetadata.samples:
        assert float(sample.sixteens_full.avgdepth['gi|219846739|ref|NR_026331.1|']) >= 15.82


def test_gdcs():
    method.run_gdcs()
    for sample in method.runmetadata.samples:
        assert os.path.isfile(sample.GDCS.baitedfastq)


def test_genesippr():
    method.genesippr()
    for sample in method.runmetadata.samples:
        assert float(sample.genesippr.avgdepth['VT1']) >= 17.38


def test_plasmids():
    method.plasmids()
    for sample in method.runmetadata.samples:
        assert float(sample.plasmidfinder.avgdepth['IncFIB(AP001918)_1__AP001918']) >= 98.66


def test_ressippr():
    method.ressippr()
    for sample in method.runmetadata.samples:
        assert float(sample.resfinder.avgdepth['sul1_1_AY224185']) >= 61.94


def test_resfinder():
    method.resfinder()
    for sample in method.runmetadata.samples:
        assert sample.resfinder_assembled.protseq


def test_prophages():
    method.prophages(cutoff=25)
    for sample in method.runmetadata.samples:
        assert sample.prophages.blastresults


def test_univec():
    method.univec()
    for sample in method.runmetadata.samples:
        assert sample.univec.blastresults


def test_virulence():
    method.virulence()
    for sample in method.runmetadata.samples:
        assert float(sample.virulence.avgdepth['stx1_3_M19437_3']) >= 16.54


def test_mlst():
    method.mlst()
    for sample in method.runmetadata.samples:
        assert sample.mlst.sequencetype == '11'


def test_serosippr():
    method.serosippr()
    for sample in method.runmetadata.samples:
        assert sample.serosippr.o_set == ['O157']


def test_vtyper():
    method.vtyper()
    for sample in method.runmetadata.samples:
        assert sample.vtyper.profile == ['vtx1a', 'vtx2b', 'vtx1d']


def test_coregenome():
    method.coregenome()
    for sample in method.runmetadata.samples:
        assert sample.coregenome.coreresults == '1/1'


def test_clear_results(variables):
    shutil.rmtree(os.path.join(variables.path, 'NC_002695'))


def test_clear_sistr(variables):
    shutil.rmtree(os.path.join(variables.path, 'NC_003198'))


def test_clear_confindr(variables):
    shutil.rmtree(os.path.join(variables.path, 'confindr'))


def test_clear_reports(variables):
    shutil.rmtree(os.path.join(variables.path, 'reports'))


def test_clear_assemblies(variables):
    shutil.rmtree(os.path.join(variables.path, 'BestAssemblies'))


def test_clear_logs(variables):
    os.remove(os.path.join(variables.path, 'logfile_err.txt'))
    os.remove(os.path.join(variables.path, 'logfile_out.txt'))
    os.remove(os.path.join(variables.path, 'portal.log'))

