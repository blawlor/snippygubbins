######################################################################
#  _   _ ____  _ _ _           
# | \ | / ___|(_) (_) ___ ___  
# |  \| \___ \| | | |/ __/ _ \ 
# | |\  |___) | | | | (_| (_) |
# |_| \_|____/|_|_|_|\___\___/ 
#
# Simplicity Tree Tools Image                             
#
# This file instructs Docker to build an image containing bioinformatic tools required by the Tree service.
#

######################################################################
FROM		nsilico/tree-base:latest
MAINTAINER 	NSilico Bioinformatics <simplicity@nsilico.com>
LABEL		description="Tree Tools" version="1.0" Vendor="NSilico Ltd."
######################################################################


USER	root
WORKDIR	/root

# Install ClonalFrameML
RUN git clone https://github.com/xavierdidelot/ClonalFrameML /tmp/ClonalFrameML
RUN cd /tmp/ClonalFrameML/src \
    && ./make.sh \
    && cp ClonalFrameML /usr/local/bin

# Gemma
ADD ./gemma.linux /usr/local/bin
RUN  chmod +x /usr/local/bin/gemma.linux
# RAxML
RUN git clone https://github.com/stamatak/standard-RAxML.git /tmp/raxml
RUN cd /tmp/raxml \
       && make -f Makefile.gcc \
       && cp raxmlHPC /usr/local/bin

ADD ./PhyML-3.1_linux64 /usr/local/bin/phyml

# Install openblas (needed by eigensoft)
RUN apt-get  update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
       libopenblas-dev \
       libgsl-dev \
       liblapacke-dev

RUN echo 'options(\n\
  repos = c(CRAN = "https://cran.r-project.org/"),\n\
  download.file.method = "libcurl",\n\
  # Detect number of physical cores\n\
  Ncpus = parallel::detectCores(logical=FALSE)\n\
)' >> /etc/R/Rprofile.site

RUN Rscript -e 'install.packages("devtools",dependencies=TRUE)' \
	&&Rscript -e 'install.packages("genoPlotR",dependencies=TRUE)'

# Eigensoft 6.1.4
ADD ./smartpca /usr/local/bin
RUN  chmod +x /usr/local/bin/smartpca

# GATK
ADD ./gatk-package-4.beta.3-48-g403ba13-SNAPSHOT-local.jar /usr/local/bin/gatk.jar


# Update the repository sources list
RUN apt-get update

# Install perl modules 
RUN apt-get install -y cpanminus \
                       vcftools \
                       make

RUN cpanm Bio::SeqIO 

ADD /snippy/binaries/noarch/snpEff /usr/local/bin/snpEff
ADD /snippy/binaries/noarch/snpEff.jar /usr/local/bin/snpEff.jar
ADD /snippy/etc/snpeff.config /root/snpeff.config
ADD /snippy/binaries/noarch/parallel /usr/local/bin/parallel
ADD /snippy/binaries/noarch/freebayes-parallel /usr/local/bin/freebayes-parallel
ADD /snippy/binaries/linux/freebayes /usr/local/bin/freebayes
ADD /snippy/binaries/noarch/fasta_generate_regions.py /usr/local/bin/fasta_generate_regions.py
ADD /snippy/binaries/linux/vcfstreamsort /usr/local/bin/vcfstreamsort
ADD /snippy/binaries/linux/vcfuniq /usr/local/bin/vcfuniq
ADD /snippy/binaries/noarch/vcffirstheader /usr/local/bin/vcffirstheader
ADD /snippy/binaries/noarch/vcf-consensus /usr/local/bin/vcf-consensus
ADD /snippy/bin/snippy-vcf_to_tab /usr/local/bin/snippy-vcf_to_tab
ADD /snippy/bin/snippy-vcf_report /usr/local/bin/snippy-vcf_report
ADD /snippy/bin/snippy-vcf_filter /usr/local/bin/snippy-vcf_filter
ADD /snippy/binaries/linux/samtools /usr/local/bin/samtools
ADD /snippy/binaries/linux/bwa /usr/local/bin/bwa
ADD /snippy/binaries/linux/tabix /usr/local/bin/tabix
ADD /snippy/binaries/linux/bgzip /usr/local/bin/bgzip

# Add Snippy
ADD /snippy/bin /root/bin
ADD /snippy/etc /root/etc

# Picard
ADD ./picard-2.10.10-SNAPSHOT-all.jar /usr/local/bin/picard.jar

RUN apt-get install -y gubbins

