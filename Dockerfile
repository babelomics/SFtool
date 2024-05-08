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
    python3-pip





ADD https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 bcftools-1-20.tar.bz2
RUN tar -xf bcftools-1-20.tar.bz2
WORKDIR "/bcftools-1.20"
RUN ./configure --prefix=/
RUN make
RUN make install
RUN rm /bcftools-1-20.tar.bz2

COPY ./docker_downloads/annovar.latest.tar.gz /docker_downloads/
WORKDIR "/docker_downloads"
RUN tar xvzf annovar.latest.tar.gz
RUN rm annovar.latest.tar.gz

ADD https://zenodo.org/records/8045374/files/hs37d5.genome.tgz?download=1 /docker_downloads/hs37d5.genome.tgz
RUN tar xvzf hs37d5.genome.tgz
RUN rm hs37d5.genome.tgz

ADD https://github.com/WGLab/InterVar/archive/refs/tags/v2.2.1.tar.gz /docker_downloads/InterVar-2.2.1.tar.gz
RUN tar -xf InterVar-2.2.1.tar.gz
WORKDIR "/docker_downloads/InterVar-2.2.1"
RUN sed -i 's/.\/convert2annovar.pl/\/docker_downloads\/annovar\/convert2annovar.pl/g' config.ini
RUN sed -i 's/.\/table_annovar.pl/\/docker_downloads\/annovar\/table_annovar.pl/g' config.ini
RUN sed -i 's/.\/annotate_variation.pl/\/docker_downloads\/annovar\/annotate_variation.pl/g' config.ini
ENV PATH "$PATH:/docker_downloads/InterVar-2.2.1/"
WORKDIR "/"

RUN python -m pip install pandas vcfpy biomart natsort pybedtools


WORKDIR "/"