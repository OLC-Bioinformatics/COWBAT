# Dockerfile for COWBAT
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	bash \
	nfs-common \
	nfs-client \
	alien \
	git \
	curl \
	libexpat1-dev \
	libxml2-dev \
	libxslt-dev \
	liblzma-dev \
	zlib1g-dev \
	libbz2-dev \
	software-properties-common \
	nano \
	xsltproc \
	python3-dev \
	libncurses5-dev \ 
        pkg-config \ 
        automake \
	libtool \
	build-essential \
	ncbi-blast+ \
	fastx-toolkit \
	python3-pip \
	autoconf && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*	

# Upgrade pip
RUN pip3 install --upgrade pip

# Add the scripts
ADD accessoryfiles/bin /accessoryfiles

# Add the databases
#ADD accessoryfile/databases /databases

# Install bcl2fastq
RUN alien -i /accessoryfiles/bcl2fastq-1.8.4-Linux-x86_64.rpm
# Remove the rpm
RUN rm /accessoryfiles/bcl2fastq-1.8.4-Linux-x86_64.rpm
# Edited Config.pm supplied with bcl2fastq to comment out sub _validateEland subroutine that was causing bcl2fastq to fail with compilation errors
COPY Config.pm /usr/local/lib/bcl2fastq-1.8.4/perl/Casava/Alignment/Config.pm

# Install cpan minus, XML:Simple and dependencies for bcl2fastq 
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm XML::Simple@2.24 --mirror-only --force

# Add bbmap files to the path
ENV PATH /accessoryfiles/bbmap:$PATH

# Install fastqc
ENV PATH /accessoryfiles/FastQC:$PATH

# Install SPAdes
ENV PATH /accessoryfiles/SPAdes/bin:$PATH

# Install qualimap
ENV PATH /accessoryfiles/qualimap:$PATH

# Install quast
ENV PATH /accessoryfiles/quast:$PATH

# Install samtools
WORKDIR /accessoryfiles/samtools
RUN make && \ 
    make prefix=/accessoryfiles/samtools install
ENV PATH /accessoryfiles/samtools/bin:$PATH

# Install jellyfish
WORKDIR /accessoryfiles/jellyfish
RUN ./configure
RUN make -j 4
RUN make install

# Install CLARK
WORKDIR /accessoryfiles/CLARK
RUN ./install.sh
ENV PATH /accessoryfiles/CLARK:$PATH
WORKDIR /

# Install Prodigal
ENV PATH /accessoryfiles/prodigal:$PATH

# Install mash
ENV PATH /accessoryfiles/mash:$PATH

# Install ePCR
ENV PATH /accessoryfiles/ePCR:$PATH

# Install bowtie2
ENV PATH /accessoryfiles/bowtie2:$PATH

# Install seqtk
RUN cd /accessoryfiles/seqtk && make
ENV PATH /accessoryfiles/seqtk:$PATH

# Install fastx_toolkit
RUN cd /accessoryfiles/libgtextutils && ./configure && make && make install
RUN cd /accessoryfiles/fastx_toolkit && ./configure && make && make install

# Install conda
RUN bash /accessoryfiles/miniconda/Miniconda3-latest-Linux-x86_64.sh -b -p /accessoryfiles/miniconda/miniconda
ENV PATH /accessoryfiles/miniconda/miniconda/bin:$PATH
RUN conda update conda
# Add Bioconda channel and other channels https://bioconda.github.io/
RUN conda config --add channels conda-forge
RUN conda config --add channels defaults
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels anaconda

# Install sistr_cmd and its dependencies
RUN conda install sistr_cmd==1.0.2

# Install pysam
RUN conda install -c bioconda pysam==0.13

# Install OLCTools
RUN pip3 install OLCTools==0.3.4

# Install the pipeline
RUN git clone https://github.com/OLC-Bioinformatics/COWBAT.git
ENV PATH /COWBAT:$PATH

# Install sipprverse
RUN pip3 install sipprverse==0.0.2

# Install python requirements
RUN cd /accessoryfiles/requirements && pip3 install -r requirements.txt

# Remove all the non-folders from the accessoryfiles folder
RUN find /accessoryfiles -name "*" -type f -exec rm -f {} \;
