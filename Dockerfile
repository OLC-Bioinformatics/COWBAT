# Dockerfile for COWBAT
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
RUN echo $PATH
	    # && rm -rf miniconda.sh \
RUN conda install -y python=3.6 && conda update conda	
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
# Add miniconda to the PATH
# ENV PATH $HOME/miniconda/bin:$PATH

# Upgrade pip
RUN pip install --upgrade pip

# Install the pipeline
WORKDIR /home/ubuntu/
ENV PATH /home/ubuntu/COWBAT:$PATH
RUN git clone https://github.com/OLC-Bioinformatics/COWBAT.git
WORKDIR /home/ubuntu/COWBAT
RUN git fetch --tags
RUN conda env create -f environment.yml

# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8

# Work-around to get MOB-suite dependencies to work
USER root
RUN ln -s /home/ubuntu/miniconda/envs/cowbat/bin/show-coords /usr/local/bin/show-coords
USER ubuntu

#CMD /bin/bash -c "source activate cowbat && assembly_pipeline.py -s /mnt/scratch/test/sequences -r /mnt/nas/assemblydatabases/0.3.4/databases"
