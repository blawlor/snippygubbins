FROM ubuntu:16.04

ENV PATH="/root/bin:${PATH}"

# Update the repository sources list
RUN apt-get update

# Install perl modules 
RUN apt-get install -y cpanminus \
                       vcftools \
                       make

RUN cpanm FindBin::libs \
          Bio::SeqIO \
          VCF

# Java need for snpEff

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

ADD /binaries/noarch/snpEff /usr/local/bin/snpEff
ADD /binaries/noarch/snpEff.jar /usr/local/bin/snpEff.jar
ADD /etc/snpeff.config /root/snpeff.config
ADD /binaries/noarch/parallel /usr/local/bin/parallel
ADD /binaries/noarch/freebayes-parallel /usr/local/bin/freebayes-parallel
ADD /binaries/linux/freebayes /usr/local/bin/freebayes
ADD /binaries/noarch/fasta_generate_regions.py /usr/local/bin/fasta_generate_regions.py
ADD /binaries/linux/vcfstreamsort /usr/local/bin/vcfstreamsort
ADD /binaries/linux/vcfuniq /usr/local/bin/vcfuniq
ADD /binaries/noarch/vcffirstheader /usr/local/bin/vcffirstheader
ADD /binaries/noarch/vcf-consensus /usr/local/bin/vcf-consensus
ADD /bin/snippy-vcf_to_tab /usr/local/bin/snippy-vcf_to_tab
ADD /bin/snippy-vcf_report /usr/local/bin/snippy-vcf_report
ADD /bin/snippy-vcf_filter /usr/local/bin/snippy-vcf_filter
ADD /binaries/linux/samtools /usr/local/bin/samtools
ADD /binaries/linux/bwa /usr/local/bin/bwa
ADD /binaries/linux/tabix /usr/local/bin/tabix
ADD /binaries/linux/bgzip /usr/local/bin/bgzip

ADD /snippy/bin /root/bin
ADD /snippy/etc /root/etc

# Install Python packages
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        python-dev \
        python-pip

