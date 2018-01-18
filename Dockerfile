# Dockerfile for COWBAT
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	wget \
	python3-pip \
	nano && \
	curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	    && bash /tmp/miniconda.sh -bfp /usr/local \
	    && rm -rf /tmp/miniconda.sh \
	    && conda install -y python=3 \
	    && conda update conda && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*	

# Upgrade pip
RUN pip3 install --upgrade pip

# Install bcl2fastq
RUN conda install -c dranew bcl2fastq

# Install bbmap 
RUN conda install -c bioconda bbmap

# Install fastqc
RUN conda install -c bioconda fastqc

# Install SPAdes
RUN conda install -c bioconda spades

# Install qualimap
RUN conda install -c bioconda qualimap

# Install samtools
RUN conda install -c bioconda samtools

# Install jellyfish
RUN conda install -c conda-forge jellyfish

# Install CLARK
RUN conda install -c bioconda clark

# Install Prodigal
RUN conda install -c biocore prodigal

# Install bowtie2 
RUN conda install -c bioconda bowtie2

# Install seqtk
RUN conda install -c bioconda seqtk

# Install fastx_toolkit
RUN conda install -c biobuilds fastx-toolkit

# Install sistr_cmd and its dependencies
RUN conda install -c bioconda sistr_cmd==1.0.2

# Install mash
RUN conda install -c bioconda mash

# Install pysam
RUN pip3 install pysam==0.13

# Install biopython 
RUN pip3 install biopython==1.70

# Install OLCTools
RUN pip3 install OLCTools

# Install sipprverse
RUN pip3 install sipprverse

# Install confindr
RUN pip3 install confindr

# Install pytest
RUN pip3 install -U pytest

# Install the pipeline
RUN git clone https://github.com/OLC-Bioinformatics/COWBAT.git
ENV PATH /COWBAT:$PATH
