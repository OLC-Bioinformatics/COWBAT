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
	ttf-dejavu \
	nano && \
	curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	    && bash /tmp/miniconda.sh -bfp /usr/local \
	    && rm -rf /tmp/miniconda.sh \
	    && conda install -y python=3 \
	    && conda update conda && \
    	rm -rf /var/lib/apt/lists/*	

# Add miniconda to the PATH
ENV PATH $HOME/miniconda/bin:$PATH

# Upgrade pip
RUN pip3 install --upgrade pip

#------

# Install the pipeline
ENV PATH /COWBAT:$PATH
RUN git clone https://github.com/OLC-Bioinformatics/COWBAT.git
WORKDIR /COWBAT
RUN conda env create
WORKDIR /

# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8

#CMD /bin/bash -c "source activate cowbat && assembly_pipeline.py /mnt/scratch/test/sequences -r /mnt/nas/assemblydatabases/0.2.1/databases"
