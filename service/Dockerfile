######################################################################
#  _   _ ____  _ _ _           
# | \ | / ___|(_) (_) ___ ___  
# |  \| \___ \| | | |/ __/ _ \ 
# | |\  |___) | | | | (_| (_) |
# |_| \_|____/|_|_|_|\___\___/ 
#
# Simplicity Tree Service Image
#
# This file instructs Docker to build an image containing the gwas service.
#

######################################################################
FROM		nsilico/tree-tools:latest
MAINTAINER 	NSilico Bioinformatics <simplicity@nsilico.com>
LABEL		description="Tree Service" version="1.0" Vendor="NSilico Ltd."
######################################################################


ADD ./main.py /root/code/main.py
ADD ./run.sh /root/run.sh
ENV NUM_CPUS 4
ENV USE_SNIPPY True
ENV USE_GUBBINS True
CMD ["./run.sh"]
