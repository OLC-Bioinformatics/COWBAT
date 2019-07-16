#Dockerfile for COWBAT

FROM ubuntu:18.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@canada.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	wget \
	python3-pip \
	ttf-dejavu \
	nano  

ENV PATH /usr/sbin:$PATH
RUN useradd -ms /bin/bash/ ubuntu
USER ubuntu

WORKDIR HOME
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/ubuntu/miniconda.sh
RUN bash /home/ubuntu/miniconda.sh -b -p /home/ubuntu/miniconda
ENV PATH /home/ubuntu/miniconda/bin:$PATH
RUN conda install -y python=3.6 && conda update conda
RUN conda config --add channels dranew	
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# Install the pipeline
RUN conda install -c olcbioinformatics cowbat

# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8

# Work-around to get MOB-suite dependencies to work
USER root
RUN ln -s /home/ubuntu/miniconda/envs/cowbat/bin/show-coords /usr/local/bin/show-coords
USER ubuntu

#CMD /bin/bash -c "assembly_pipeline.py -s /mnt/scratch/test/sequences -r /mnt/nas/assemblydatabases/0.5.0.0"
