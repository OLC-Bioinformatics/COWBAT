SPAdesPipeline
==============
#Introduction

This pipeline is designed to be able to assemble and type raw fastq data generated using an Illumina MiSeq.

It has been designed to run using either the archived files obtained from BaseSpace (something similar to analysis_14348334_fastq.zip),
alternatively, it is possible to use the individual fastq files taken directly from the MiSeq. This is probably better, as,
in addition to the fastq files, three other files from the MiSeq are required for the pipeline to function:

1. GenerateFASTQRunStatistics.xml
2. RunInfo.xml
3. SampleSheet.csv

These files are located within the appropriate subfolder (e.g. 140922_M02466_0030_000000000-AARWU - this folder consists
of the date (140922), the MiSeq designation (M02466), and the flowcell number (000000000-AARWU)) of the MiSeqOutput
directory in the MiSeq onboard computer. The fastq files are located in the ../MiSeqOutput/140922_M02466_0030_000000000-AARWU/Data/Intensities/BaseCalls
folder.

Copy all necessary files to a properly named folder in an easy to remember location (e.g. ../Sequencing/140922).

#Contents
This pipeline includes a main script (SPAdesPipeline), and several helper modules located in helper scripts, including
* Pulling metadata from sequencing run reports/files
 * runMetadataOptater
* Moving and/or extracting archived files
 * fileExtractionProcessing
* Performing quake error corrections on the reads
 * quakeR
* Running SPAdes assembler
 * spadesGoUpper
* Performing typing using rMLST to determine best reference genome
 * rMLST_typer
* Determining assembly quality metrics using quast
 * quastR
* Estimating the size of the library fragments
 * lse
* Creating a JSON report of all collected metadata for each sequenced strain
 * reportR

#Use
Run SPAdesPipeline.py from the console.

#Requirements
* Docker
Everything else should be contained within the docker container, and is ready to run.

# Outputs
This pipeline generates multiple outputs.

1. Assembled contigs - these are collected in the 'BestAssemblies' folder
2. JSON reports - these are located in the 'reports' folder
3. A summary of rMLST alleles - this is located in the 'reports' folder

Additionally, within the individual strain subfolders, a .pdf output of plotted insert sizes is included in the 'insertSizes' folder.
Detailed reports can be found in the 'quast_results' folder, and the reference genome file is located in 'referenceGenome'

