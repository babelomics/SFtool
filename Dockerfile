FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Madrid

RUN apt-get update \
    && apt install -y \
    build-essential \
    checkinstall \
    curl \
    cmake \
    git \
    help2man \
    lsb-release \
    rpm \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3.12 \
    python3-pip \
    locales \
    wget \
    zip \
    gzip \
    samtools

RUN sed -i '/es_ES.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG es_ES.UTF-8
ENV LANGUAGE es_ES.UTF-8
ENV LC_ALL es_ES.UTF-8

RUN mkdir -p /docker_files
COPY ./docker_files/config_docker.json /docker_files

RUN mkdir -p /docker_directories/ref_genomes/38
ADD http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz /docker_directories/ref_genomes/38/hg38.fa.gz
WORKDIR "/docker_directories/ref_genomes/38/"
RUN gunzip hg38.fa.gz
RUN samtools faidx hg38.fa


RUN mkdir -p /docker_dependencies
WORKDIR "/docker_dependencies"
ADD https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 bcftools-1.21.tar.bz2
RUN tar -xf bcftools-1.21.tar.bz2
WORKDIR "/docker_dependencies/bcftools-1.21"
RUN ./configure --prefix=/
RUN make
RUN make install
RUN rm /docker_dependencies/bcftools-1.21.tar.bz2

WORKDIR "/docker_dependencies"
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz bedtools-2.31.1.tar.gz
RUN tar -xf bedtools-2.31.1.tar.gz
WORKDIR "/docker_dependencies/bedtools2"
RUN make
ENV PATH "$PATH:/docker_dependencies/bedtools2/bin/"
WORKDIR "/docker_dependencies"
RUN rm bedtools-2.31.1.tar.gz

WORKDIR "/docker_dependencies"
ADD https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 htslib-1.21.tar.bz2
RUN tar -xf htslib-1.21.tar.bz2
WORKDIR "/docker_dependencies/htslib-1.21"
RUN ./configure --prefix=/
RUN make
RUN make install
RUN rm /docker_dependencies/htslib-1.21.tar.bz2

WORKDIR "/docker_dependencies"
ADD https://github.com/PharmGKB/PharmCAT/releases/download/v2.15.5/pharmcat-pipeline-2.15.5.tar.gz pharmcat-pipeline-2.15.5.tar.gz
RUN tar -xf pharmcat-pipeline-2.15.5.tar.gz
RUN rm pharmcat-pipeline-2.15.5.tar.gz

WORKDIR "/docker_dependencies"
ADD https://github.com/pstawinski/genebe-cli/releases/download/v0.1.0-a.6/GeneBeClient-0.1.0-a.6.jar GeneBeClient-0.1.0-a.6.jar
RUN ln -s GeneBeClient-0.1.0-a.6.jar GeneBeClient.jar


RUN mkdir -p /docker_directories/ref_genomes/37
ADD https://zenodo.org/records/8045374/files/hs37d5.genome.tgz?download=1 /docker_directories/ref_genomes/37/hs37d5.genome.tgz
WORKDIR "/docker_directories/ref_genomes/37/"
RUN tar xvzf hs37d5.genome.tgz
RUN rm hs37d5.genome.tgz



RUN pip3 install pandas vcfpy requests natsort pybedtools openpyxl --break-system-packages

RUN mkdir -p /docker_directories/categories/
COPY categories /docker_directories/categories

ADD https://zenodo.org/records/15068802/files/clinvar_database_GRCh37_20250307.txt?download=1 /docker_dependencies/clinvar/clinvar_database_GRCh37_20250307.txt
WORKDIR "/docker_dependencies/clinvar/"
RUN chmod 755 clinvar_database_GRCh37_20240528.txt

ADD https://zenodo.org/records/15068811/files/clinvar_database_GRCh38_20250307.txt?download=1 /docker_dependencies/clinvar/clinvar_database_GRCh38_20250307.txt
WORKDIR "/docker_dependencies/clinvar/"
RUN chmod 755 clinvar_database_GRCh38_20250307.txt


RUN mkdir -p /release_build/

WORKDIR "."
COPY modules /release_build/modules
COPY SFtool.py /release_build/SFtool.py
COPY docker_files/SFtool_singularity /release_build/SFtool
ENV PATH "$PATH:/release_build/"

WORKDIR "/"