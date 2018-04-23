#####################################################################e
#  _   _ ____  _ _ _           
# | \ | / ___|(_) (_) ___ ___  
# |  \| \___ \| | | |/ __/ _ \ 
# | |\  |___) | | | | (_| (_) |
# |_| \_|____/|_|_|_|\___\___/ 
#
# Simplicity Base Image for Tree
#
# This file instructs Docker to build a base image to be used by the Profiler Tools image.
# 
######################################################################
FROM		ubuntu:16.04
MAINTAINER 	NSilico Bioinformatics <simplicity@nsilico.com>

######################################################################

# Select an Irish mirror
#RUN 	sed 's@archive.ubuntu.com@ftp.heanet.ie/pub/@' -i /etc/apt/sources.list

# Install system packages

RUN	apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
	aria2 \
	ca-certificates \
	cron-apt \
	g++ \
	libcurl4-openssl-dev \
	make \
	unzip \
	vim \
	wget \
        git \
	xorg-dev \
        software-properties-common \
        apt-transport-https

RUN     DEBIAN_FRONTEND=noninteractive apt-get install -y libncurses5-dev

RUN     gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN     gpg -a --export E084DAB9 | apt-key add -

RUN     echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
 
RUN     apt-get update

# Install R evnvironment
RUN     apt-get install -y r-base-core r-recommended r-base r-base-dev

# This is in accordance to : https://www.digitalocean.com/community/tutorials/how-to-install-java-with-apt-get-on-ubuntu-16-04
RUN apt-get update && \
	apt-get install -y openjdk-8-jdk && \
	apt-get install -y ant && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;
	
# Fix certificate issues, found as of 
# https://bugs.launchpad.net/ubuntu/+source/ca-certificates-java/+bug/983302
RUN apt-get update && \
	apt-get install -y ca-certificates-java && \
	apt-get clean && \
	update-ca-certificates -f && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer;

# Setup JAVA_HOME, this is useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

# Install Python packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
	python-dev \
	python-pip

RUN	pip install --upgrade distribute && \
	pip install biopython pyvcf


