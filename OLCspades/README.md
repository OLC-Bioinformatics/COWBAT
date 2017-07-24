OLC SPAdes Assembly Pipeline
===========================
#Introduction

This pipeline is designed to be able to assemble and type raw fastq data generated using an Illumina MiSeq.

It has been designed to run using either the archived files obtained from BaseSpace (something similar to analysis_14348334_fastq.zip),
alternatively, it is possible to use the individual fastq files taken directly from the MiSeq. This is probably better, as,
in addition to the fastq files, three other files from the MiSeq are required for the pipeline to function with full metadata retrieval:

1. GenerateFASTQRunStatistics.xml
2. RunInfo.xml
3. SampleSheet.csv

These files are located within the appropriate subfolder (e.g. 140922_M02466_0030_000000000-AARWU - this folder consists
of the date (140922), the MiSeq designation (M02466), and the flowcell number (000000000-AARWU)) of the MiSeqOutput
directory in the onboard MiSeq computer. The fastq files are located in the ../MiSeqOutput/140922_M02466_0030_000000000-AARWU/Data/Intensities/BaseCalls
folder.

Copy all necessary files to a properly named folder in an easy to remember location (e.g. ../Sequencing/140922).

#Requirements
* bbduk.sh and spades.py in $PATH
* Docker
Everything else should be contained within the docker container, and is ready to run.

#Contents
This pipeline includes a main script (OLCspades.py), and several helper modules located in helper scripts.

#Use
Run OLCspades.py from the console.


```
usage: OLCspades.py [-h] [-v] [-n NUMREADS] [-t THREADS] [-o] [-F]
                    [-d DESTINATIONFASTQ] [-m MISEQPATH] [-f MISEQFOLDER]
                    [-r1 READLENGTHFORWARD] [-r2 READLENGTHREVERSE]
                    [-r REFERENCEFILEPATH] [-k KMERRANGE]
                    [-c CUSTOMSAMPLESHEET] [-b] [-p] [-u]
                    path

Assemble genomes from Illumina fastq files

positional arguments:
  path                  Specify path

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n NUMREADS, --numreads NUMREADS
                        Specify the number of reads. Paired-reads: 2,
                        unpaired-reads: 1. Default is paired-end
  -t THREADS, --threads THREADS
                        Number of threads. Default is the number of cores in
                        the system
  -o, --offhours        Optionally run the off-hours module that will search
                        for MiSeq runs in progress, wait until the run is
                        complete, and assemble the run
  -F, --fastqcreation   Optionally run the fastq creation module that will
                        search for MiSeq runs in progress, run bcl2fastq to
                        create fastq files, and assemble the run
  -d DESTINATIONFASTQ, --destinationfastq DESTINATIONFASTQ
                        Optional folder path to store .fastq files created
                        using the fastqCreation module. Defaults to
                        path/miseqfolder
  -m MISEQPATH, --miseqpath MISEQPATH
                        Path of the folder containing MiSeq run data folder
  -f MISEQFOLDER, --miseqfolder MISEQFOLDER
                        Name of the folder containing MiSeq run data
  -r1 READLENGTHFORWARD, --readlengthforward READLENGTHFORWARD
                        Length of forward reads to use. Can specify "full" to
                        take the full length of forward reads specified on the
                        SampleSheet. Defaults to "full"
  -r2 READLENGTHREVERSE, --readlengthreverse READLENGTHREVERSE
                        Length of reverse reads to use. Can specify "full" to
                        take the full length of reverse reads specified on the
                        SampleSheet. Defaults to "full"
  -r REFERENCEFILEPATH, --referencefilepath REFERENCEFILEPATH
                        Provide the location of the folder containing the
                        pipeline accessory files (reference genomes, MLST
                        data, etc.
  -k KMERRANGE, --kmerrange KMERRANGE
                        The range of kmers used in SPAdes assembly. Default is
                        21,33,55,77,99,127
  -c CUSTOMSAMPLESHEET, --customsamplesheet CUSTOMSAMPLESHEET
                        Path of folder containing a custom sample sheet and
                        name of sample sheet file e.g.
                        /home/name/folder/BackupSampleSheet.csv. Note that
                        this sheet must still have the same format of Illumina
                        SampleSheet.csv files
  -b, --basicassembly   Performs a basic de novo assembly, and does not
                        collect run metadata
  -p, --preprocess      Perform quality trimming and error correction only. Do
                        not assemble the trimmed + corrected reads
  -u, --updatedatabases
                        Optionally update (r)MLST databases
```


# Outputs
This pipeline generates multiple outputs.

1. Assembled contigs - these are collected in the 'BestAssemblies' folder
2. JSON reports - these are located in the 'reports' folder
3. A summary of rMLST alleles - this is located in the 'reports' folder


## Things to add to this README:
* bbduk.sh and spades.py need to be in your $PATH
* SPAdes needs to be at a version that supports python3 - using 3.10.1 seems to work.
* BBduk/the bbmap package needs to be upgraded to version 37.23
* Quast needs to be upgraded to version 4.5 in order to support python3

