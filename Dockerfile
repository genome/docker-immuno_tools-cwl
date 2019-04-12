FROM fred2/optitype:latest
#unfortunately the most up-to-date version is only tagged latest, no release tag; however, it hasn't been updated in a year
MAINTAINER John Garza <johnegarza@wustl.edu>

LABEL \
    description="Image containing optitype and other helper tools for the optitype HLA-typing immuno subworkflow"

#switch from inherited user to root for necessary installation permissions
USER root

#needed for samtools
RUN apt-get update -y && apt-get install -y libncurses5-dev

#htslib, samtools, sambamba, bwa install below all pulled from docker-cle dockerfile

##############
#HTSlib 1.3.2#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /tmp/htslib-1.3.2
RUN ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

################
#Samtools 1.3.1#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar --bzip2 -xf samtools-1.3.1.tar.bz2

WORKDIR /tmp/samtools-1.3.1
RUN ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install

WORKDIR /
RUN rm -rf /tmp/samtools-1.3.1

#################
#Sambamba v0.6.4#
#################

RUN mkdir /opt/sambamba/ \
    && wget https://github.com/lomereiter/sambamba/releases/download/v0.6.4/sambamba_v0.6.4_linux.tar.bz2 \
    && tar --extract --bzip2 --directory=/opt/sambamba --file=sambamba_v0.6.4_linux.tar.bz2 \
    && ln -s /opt/sambamba/sambamba_v0.6.4 /usr/bin/sambamba

############
#BWA 0.7.15#
############

ENV BWA_VERSION 0.7.15

RUN cd /tmp/ \
    && wget -q http://downloads.sourceforge.net/project/bio-bwa/bwa-${BWA_VERSION}.tar.bz2 && tar xvf bwa-${BWA_VERSION}.tar.bz2 \
    && cd /tmp/bwa-${BWA_VERSION} \
    && sed -i 's/CFLAGS=\\t\\t-g -Wall -Wno-unused-function -O2/CFLAGS=-g -Wall -Wno-unused-function -O2 -static/' Makefile \
    && make \
    && cp /tmp/bwa-${BWA_VERSION}/bwa /usr/local/bin \
    && rm -rf /tmp/bwa-${BWA_VERSION}

#################
#bedtools 2.28.0#
#################
WORKDIR /usr/bin/
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools
RUN chmod +x bedtools

COPY optitype_script.sh /usr/bin/optitype_script.sh

#include the optitype hla reference and pre-calculated bwa index outputs
RUN mkdir -p /ref_data/optitype_ref
RUN mkdir -p /ref_data/ebi_ref
COPY hla_reference_dna.fasta /ref_data/optitype_ref/hla_reference_dna.fasta
COPY hla_reference_dna.fasta.amb /ref_data/optitype_ref/hla_reference_dna.fasta.amb
COPY hla_reference_dna.fasta.ann /ref_data/optitype_ref/hla_reference_dna.fasta.ann
COPY hla_reference_dna.fasta.bwt /ref_data/optitype_ref/hla_reference_dna.fasta.bwt
COPY hla_reference_dna.fasta.pac /ref_data/optitype_ref/hla_reference_dna.fasta.pac
COPY hla_reference_dna.fasta.sa /ref_data/optitype_ref/hla_reference_dna.fasta.sa
COPY hla_nuc.fasta /ref_data/ebi_ref/hla_nuc.fasta
COPY hla_nuc.fasta.amb /ref_data/ebi_ref/hla_nuc.fasta.amb
COPY hla_nuc.fasta.ann /ref_data/ebi_ref/hla_nuc.fasta.ann
COPY hla_nuc.fasta.bwt /ref_data/ebi_ref/hla_nuc.fasta.bwt
COPY hla_nuc.fasta.pac /ref_data/ebi_ref/hla_nuc.fasta.pac
COPY hla_nuc.fasta.sa /ref_data/ebi_ref/hla_nuc.fasta.sa

#clear inherited entrypoint
ENTRYPOINT []
CMD []
