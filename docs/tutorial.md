## Tutorial

### Basic settings

The pipeline is designed to run with a minimum of two supplied parameters:

    * path to FASTQ sequence data (-s)
    * path to reference database (-r)

The following command will run the pipeline on the supplied sequences with default parameters
    
```
assembly_pipeline.py -s /path/to/sequences -r /path/to/database

```

### Optional parameters

There are a number of optional parameters than can be supplied to the assembly_pipeline.py script

```
  -h, --help            
                        show help message and exit
  -v, --version         
                        show program's version number and exit
  
  -n, --numreads NUMREADS
                        Specify the number of reads. Paired-reads: 2,
                        unpaired-reads: 1. Default is paired-end
  -t, --threads THREADS
                        Number of threads. Default is the number of cores in
                        the system
  -c, --customsamplesheet CUSTOMSAMPLESHEET
                        Path of folder containing a custom sample sheet and
                        name of sample sheet file e.g.
                        /home/name/folder/BackupSampleSheet.csv. Note that
                        this sheet must still have the same format of Illumina
                        SampleSheet.csv files
  -b, --basicassembly   
                        Performs a basic de novo assembly, and does not
                        collect run metadata
  -p, --preprocess      
                        Performs quality trimming and error correction only. Do
                        not assemble the trimmed + corrected reads
```