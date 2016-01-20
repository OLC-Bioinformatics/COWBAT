# Dockerfile for OLCspades genome assembly pipeline
FROM ubuntu:14.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

COPY sources.list /etc/apt/sources.list

# Install various required softwares
RUN apt-get update -y -qq && apt-get install -y --force-yes \
	bash \
	nfs-common \
	nfs-client \
	alien \
	git \
	curl \
	libexpat1-dev \
	libxml2-dev \
	libxslt-dev \
	zlib1g-dev \
	libbz2-dev \
	software-properties-common \
	nano \
	xsltproc
RUN echo 'hello'

# Install bcl2fastq
ADD accessoryfiles /accessoryfiles
RUN alien -i /accessoryfiles/bcl2fastq-1.8.4-Linux-x86_64.rpm
# Remove the rpm
RUN rm /accessoryfiles/bcl2fastq-1.8.4-Linux-x86_64.rpm
# Edited Config.pm supplied with bcl2fastq to comment out sub _validateEland subroutine that was causing bcl2fastq to fail with compilation errors
COPY Config.pm /usr/local/lib/bcl2fastq-1.8.4/perl/Casava/Alignment/Config.pm

# Install XML:Simple and dependencies for bcl2fastq 
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm --mirror http://mirror.csclub.uwaterloo.ca/CPAN/ XML::Simple --mirror-only --force

# Install bbmap and bbduk
RUN echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | sudo /usr/bin/debconf-set-selections
RUN add-apt-repository -y ppa:webupd8team/java
RUN apt-get update -qq -y && apt-get install -y --force-yes \
	oracle-java7-installer \
	oracle-java7-set-default  && \
    	rm -rf /var/cache/oracle-jdk7-installer  && \
    	apt-get clean  && \
    	rm -rf /var/lib/apt/lists/*

# Add bbmap files to the path
ENV PATH /accessoryfiles/bbmap:$PATH

# Install fastqc
ENV PATH /accessoryfiles/FastQC:$PATH

# Install SPAdes
ENV PATH /accessoryfiles/spades:$PATH

# Mount NAS:
# must run container in privileged mode (--privileged)
RUN locale-gen en_US en_US.UTF-8
RUN dpkg-reconfigure locales

COPY mount_nfs.sh /root/mount_nfs.sh

# run this script in your cmd or entrypoint script to mount your nfs mounts
RUN chmod +x /root/mount_nfs.sh
ENTRYPOINT ["/root/mount_nfs.sh"]

# Useful commands
#docker build -t remotepythondocker .
#docker run -e NFS_MOUNT=192.168.1.18:/mnt/zvolume1 --privileged -it -v /home/blais/PycharmProjects/SPAdesPipeline:/spades -v /media/miseq/:/media/miseq -v /home/blais/Downloads/accessoryfiles:/accessoryfiles --name pythondocker remotepythondocker
#docker rm pythondocker

# python /spades/OLCspades/OLCspades.py -m /media/miseq/MiSeqOutput -F /mnt/zvolume1/akoziol/Pipeline_development/OLCspadesV2 -f 151218_M02466_0126_000000000-AKF4P -r1 25 -r2 25 -r /mnt/zvolume1/akoziol/WGS_Spades/spades_pipeline/SPAdesPipelineFiles/
