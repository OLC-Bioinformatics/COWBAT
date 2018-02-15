## Tutorial

### Basic settings

The pipeline is designed to run with a minimum of two supplied parameters:

    * path to FASTQ sequence data
    * path to reference database (-r)

The following command will run the pipeline on the supplied sequences with default parameters
    
```
assembly_pipeline.py /path/to/sequences -r /path/to/database

```

### Optional parameters

There are a number of optional parameters than can be supplied to the assembly_pipeline.py script

```
usage: assembly_pipeline.py [-h] [-v] [-n NUMREADS] [-t THREADS]
                            [-k KMERRANGE] [-c CUSTOMSAMPLESHEET] [-b] [-p]

Assemble genomes from Illumina fastq files

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n NUMREADS, --numreads NUMREADS
                        Specify the number of reads. Paired-reads: 2,
                        unpaired-reads: 1. Default is paired-end
  -t THREADS, --threads THREADS
                        Number of threads. Default is the number of cores in
                        the system
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
```