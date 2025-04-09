################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="RNAsplicing-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for RNAsplicing-nf"
LABEL about.home="http://github.com/IARCbioinfo/RNAsplicing-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/RNAsplicing-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/RNAsplicing-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@iarc.who.int**>

################## INSTALLATION ######################
COPY environment.yml /
#RUN apt-get update && apt-get install -y procps && apt-get clean -y
#RUN conda config --set channel_priority strict
#RUN conda env create -n rnaseq-nf -f /environment.yml && conda clean -a
#ENV PATH /opt/conda/envs/rnaseq-nf/bin:$PATH
RUN conda env update -n root -f /environment.yml && conda clean -a
RUN cd /
RUN git clone https://github.com/comprna/SUPPA