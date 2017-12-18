#!/usr/bin/env python 3
from Bio.Sequencing.Applications import SamtoolsFaidxCommandline, SamtoolsIndexCommandline, \
    SamtoolsSortCommandline, SamtoolsViewCommandline
from sipprCommon.bowtie import Bowtie2CommandLine, Bowtie2BuildCommandLine
from Bio.Blast.Applications import NcbiblastnCommandline
from argparse import ArgumentParser
from subprocess import call
import multiprocessing
from time import time
import logging
import pytest
import shutil
import sys
import os

testpath = os.path.abspath(os.path.dirname(__file__))
scriptpath = os.path.join(testpath, '..')
sys.path.append(scriptpath)
from method import Method

__author__ = 'adamkoziol'

# Set global variables
logging.basicConfig(level=logging.DEBUG)
runmetadata = list()


@pytest.fixture()
def variables():
    v = ArgumentParser()
    v.path = os.path.join(testpath, 'testdata', 'results')
    v.targetpath = os.path.join(testpath, 'testdata', 'targets')
    v.miseqpath = os.path.join(testpath, 'testdata')
    v.miseqfolder = 'flowcell'
    v.readlengthforward = '1'
    v.readlengthreverse = '0'
    v.customsamplesheet = os.path.join(v.miseqpath, v.miseqfolder, 'SampleSheet.csv')
    v.copy = True
    v.debug = True
    return v


@pytest.fixture()
def method_init(variables):
    method = Method(variables, '', time(), scriptpath)
    return method


method = method_init(variables())


def test_bcl2fastq(variables):
    method.createobjects()
    assert os.path.isfile(os.path.join(variables.path, variables.miseqfolder, '1_0',
                                       'Undetermined_S0_L001_R1_001.fastq.gz'))

def metadata_update(analysistype):
    """

    :param analysistype:
    :return:
    """
    method.sequencepath = os.path.join(testpath, 'testdata', 'sequences', analysistype)
    method.reportpath = os.path.join(testpath, 'testdata', 'reports')
    for sample in method.runmetadata.samples:
        sample.name = 'unit_test'
        sample.general.outputdirectory = method.sequencepath
        sample.run.outputdirectory = method.sequencepath
        sample.general.fastqfiles = [os.path.join(method.sequencepath, 'reads.fastq.gz')]
        sample.general.trimmedcorrectedfastqfiles = sample.general.fastqfiles
        sample.general.logout = os.path.join(method.sequencepath, 'logout')
        sample.general.logerr = os.path.join(method.sequencepath, 'logerr')


def test_fastq_bait(variables):
    outfile = os.path.join(variables.path, 'bait', 'baited.fastq')
    targetpath = os.path.join(variables.targetpath, 'bait')
    baitcall = 'bbduk.sh ref={ref} in={input} threads={cpus} outm={out}'.format(
        ref=os.path.join(targetpath, 'combinedtargets.fasta'),
        input=os.path.join(targetpath, 'genesippr.fastq.gz'),
        cpus=multiprocessing.cpu_count(),
        out=os.path.join(outfile)
    )
    call(baitcall, shell=True)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_reverse_bait(variables):
    outfile = os.path.join(variables.path, 'reverse_bait', 'baited_targets.fasta')
    targetpath = os.path.join(variables.targetpath, 'bait')
    baitcall = 'bbduk.sh ref={ref} in={input} threads={cpus} outm={out}'.format(
        ref=os.path.join(targetpath, 'genesippr.fastq.gz'),
        input=os.path.join(targetpath, 'combinedtargets.fasta'),
        cpus=multiprocessing.cpu_count(),
        out=os.path.join(outfile)
    )
    call(baitcall, shell=True)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_bowtie2_build(variables):
    # Use bowtie2 wrapper to create index the target file
    targetpath = os.path.join(variables.targetpath, 'bait')
    bowtie2build = Bowtie2BuildCommandLine(reference=os.path.join(targetpath, 'baitedtargets.fa'),
                                           bt2=os.path.join(targetpath, 'baitedtargets'))

    bowtie2build()
    size = os.stat(os.path.join(targetpath, 'baitedtargets.1.bt2'))
    assert size.st_size > 0


def test_bowtie2_align(variables):
    outpath = os.path.join(variables.path, 'bait')
    outfile = os.path.join(outpath, 'map_test_sorted.bam')
    targetpath = os.path.join(variables.targetpath, 'bait')
    # Use samtools wrapper to set up the bam sorting command
    samsort = SamtoolsSortCommandline(input=outfile,
                                      o=True,
                                      out_prefix="-")
    samtools = [
        # When bowtie2 maps reads to all possible locations rather than choosing a 'best' placement, the
        # SAM header for that read is set to 'secondary alignment', or 256. Please see:
        # http://davetang.org/muse/2014/03/06/understanding-bam-flags/ The script below reads in the stdin
        # and subtracts 256 from headers which include 256
        'python3 {}'.format(scriptpath),
        # Use samtools wrapper to set up the samtools view
        SamtoolsViewCommandline(b=True,
                                S=True,
                                h=True,
                                input_file="-"),
        samsort]
    # Add custom parameters to a dictionary to be used in the bowtie2 alignment wrapper
    indict = {'--very-sensitive-local': True,
              '-U': os.path.join(targetpath, 'genesippr.fastq.gz'),
              '-a': True,
              '--threads': multiprocessing.cpu_count(),
              '--local': True}
    # Create the bowtie2 reference mapping command
    bowtie2align = Bowtie2CommandLine(bt2=os.path.join(targetpath, 'baitedtargets'),
                                      threads=multiprocessing.cpu_count(),
                                      samtools=samtools,
                                      **indict)
    bowtie2align(cwd=outpath)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_index_target(variables):
    targetpath = os.path.join(variables.targetpath, 'bait')
    target_index = SamtoolsFaidxCommandline(reference=os.path.join(targetpath, 'baitedtargets.fa'))
    target_index()
    size = os.stat(os.path.join(targetpath, 'baitedtargets.fa.fai'))
    assert size.st_size > 0


def test_index_bam(variables):
    targetpath = os.path.join(variables.targetpath, 'bait')
    bam_index = SamtoolsIndexCommandline(input=os.path.join(targetpath, 'genesippr_sorted.bam'))
    bam_index()
    size = os.stat(os.path.join(targetpath, 'genesippr_sorted.bam.bai'))
    assert size.st_size > 0


def test_subsample(variables):
    targetpath = os.path.join(variables.targetpath, 'blast')
    outpath = os.path.join(variables.path, 'blast')
    os.mkdir(outpath)
    outfile = os.path.join(outpath, 'subsampled_reads.fastq.gz')
    cmd = 'reformat.sh in={input} out={output} samplebasestarget=100000'.format(
        input=os.path.join(targetpath, 'reads.fastq.gz'),
        output=os.path.join(outfile))
    call(cmd, shell=True)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_downsample(variables):
    outpath = os.path.join(variables.path, 'blast')
    outfile = os.path.join(outpath, 'subsampled_reads.fastq')
    cmd = 'seqtk sample {input} 1000 > {output}' .format(
        input=os.path.join(outpath, 'subsampled_reads.fastq.gz'),
        output=outfile)
    call(cmd, shell=True)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_fastq_to_fasta(variables):
    outfile = os.path.join(variables.path, 'blast', 'subsampled_reads.fasta')
    cmd = 'fastq_to_fasta -i {input} -o {output}'.format(
        input=os.path.join(os.path.join(variables.path, 'blast', 'subsampled_reads.fastq')),
        output=outfile)
    call(cmd, shell=True)
    size = os.stat(outfile)
    assert size.st_size > 0


def test_make_blastdb(variables):
    targetpath = os.path.join(variables.targetpath, 'blast')
    command = 'makeblastdb -in {targets} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {output}'.format(
        targets=os.path.join(targetpath, 'baitedtargets.fa'),
        output=os.path.join(targetpath, 'baitedtargets'))
    call(command, shell=True)
    outfile = os.path.join(targetpath, 'baitedtargets.nsi')
    size = os.stat(outfile)
    assert size.st_size > 0


def test_blast(variables):
    targetpath = os.path.join(variables.targetpath, 'blast')
    outpath = os.path.join(variables.path, 'blast')
    outfile = os.path.join(outpath, 'blast_results.csv')
    # Use the NCBI BLASTn command line wrapper module from BioPython to set the parameters of the search
    blastn = NcbiblastnCommandline(query=os.path.join(outpath, 'subsampled_reads.fasta'),
                                   db=os.path.join(targetpath, 'baitedtargets'),
                                   max_target_seqs=1,
                                   num_threads=multiprocessing.cpu_count(),
                                   outfmt="'6 qseqid sseqid positive mismatch gaps "
                                          "evalue bitscore slen length qstart qend qseq sstart send sseq'",
                                   out=outfile)
    blastn()
    size = os.stat(outfile)
    assert size.st_size > 0


def clean_folder(analysistype):
    """

    :param analysistype:
    """
    shutil.rmtree(os.path.join(method.sequencepath, analysistype))
    os.remove(os.path.join(method.sequencepath, 'logout'))
    os.remove(os.path.join(method.sequencepath, 'logerr'))
    os.remove(os.path.join(method.sequencepath, 'unit_test_metadata.json'))


def test_genesippr():
    analysistype = 'genesippr'
    metadata_update(analysistype)
    method.run_genesippr()
    outfile = os.path.join(method.reportpath, '{}.csv'.format(analysistype))
    size = os.stat(outfile)
    clean_folder(analysistype)
    assert size.st_size > 0


def test_sixteens():
    analysistype = 'sixteens_full'
    metadata_update(analysistype)
    method.run_sixteens()
    outfile = os.path.join(method.reportpath, '{}.csv'.format(analysistype))
    size = os.stat(outfile)
    clean_folder(analysistype)
    assert size.st_size > 0


def test_gdcs():
    analysistype = 'GDCS'
    metadata_update(analysistype)
    method.run_gdcs()
    outfile = os.path.join(method.reportpath, '{}.csv'.format(analysistype))
    size = os.stat(outfile)
    clean_folder(analysistype)
    assert size.st_size > 0

# def test_serosippr():
#     metadata_update('serosippr')
#     method.run_serosippr()


def test_clear_results(variables):
    shutil.rmtree(variables.path)


def test_clear_reports():
    shutil.rmtree(os.path.join(testpath, 'testdata', 'reports'))


def test_clear_targets(variables):
    targetpath = os.path.join(variables.targetpath, 'bait')
    os.remove(os.path.join(targetpath, 'baitedtargets.1.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.2.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.3.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.4.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.rev.1.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.rev.2.bt2'))
    os.remove(os.path.join(targetpath, 'baitedtargets.fa.fai'))
    os.remove(os.path.join(targetpath, 'genesippr_sorted.bam.bai'))


def test_clear_blast(variables):
    targetpath = os.path.join(variables.targetpath, 'blast')
    os.remove(os.path.join(targetpath, 'baitedtargets.nsq'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nsi'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nsd'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nog'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nni'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nnd'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nin'))
    os.remove(os.path.join(targetpath, 'baitedtargets.nhr'))
