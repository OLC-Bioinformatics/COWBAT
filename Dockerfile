# syntax=docker/dockerfile:1.3

# Create an image to install the COWBAT pipeline
FROM ubuntu:22.04 AS setup

# Install packages
RUN apt update -y -qq && apt install -y wget

# Change dir
WORKDIR /opt

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh && \
	bash /opt/miniconda.sh -b -p /opt/miniconda && \
	rm miniconda.sh

# Add miniconda to PATH
ENV PATH /opt/miniconda/bin:$PATH

# Add conda channels
RUN conda config --add channels Freenome && \
	conda config --add channels bioconda && \
	conda config --add channels conda-forge

# Install mamba
RUN conda install mamba -n base -c conda-forge

# Install the COWBAT pipeline with mamba
RUN mamba create -n cowbat -c olcbioinformatics -y cowbat=0.5.0.23=py_3


# Create a lightweight image with only the cowbat environment created above
FROM ubuntu:22.04

MAINTAINER Adam Koziol <adam.koziol@inspection.gc.ca>

# Copy the environment
COPY --from=setup /opt/miniconda/envs/cowbat /opt/miniconda/envs/cowbat

# Add binaries to PATH
ENV PATH /opt/miniconda/envs/cowbat/bin:$PATH

# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8

# Work-around to get MOB-suite dependencies to work
RUN ln -s /opt/miniconda/envs/cowbat/bin/show-coords /usr/local/bin/show-coords

# Set the path to the python executable
ENV PYTHONPATH=/opt/miniconda/envs/cowbat/bin/python

