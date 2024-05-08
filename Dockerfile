FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Madrid

RUN apt-get update \
    && apt install -y \
    build-essential \
    curl \
    cmake \
    git \
    help2man \
    lsb-release \
    python3 \
    python3-pip \
    rpm \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev


ADD https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 bcftools-1-20.tar.bz2
RUN tar -xf bcftools-1-20.tar.bz2
WORKDIR "/bcftools-1.20"
RUN ./configure --prefix=/
RUN make
RUN make install
RUN rm /bcftools-1-20.tar.bz2

COPY ./docker_dependencies/annovar.latest.tar.gz /docker_dependencies/
WORKDIR "/docker_dependencies"
RUN tar xvzf annovar.latest.tar.gz
RUN rm annovar.latest.tar.gz

ADD https://zenodo.org/records/8045374/files/hs37d5.genome.tgz?download=1 /docker_dependencies/hs37d5.genome.tgz
RUN tar xvzf hs37d5.genome.tgz
RUN rm hs37d5.genome.tgz

ADD https://github.com/WGLab/InterVar/archive/refs/tags/v2.2.1.tar.gz /docker_dependencies/InterVar-2.2.1.tar.gz
RUN tar -xf InterVar-2.2.1.tar.gz
WORKDIR "/docker_dependencies/InterVar-2.2.1"
RUN sed -i 's/.\/convert2annovar.pl/\/docker_dependencies\/annovar\/convert2annovar.pl/g' config.ini
RUN sed -i 's/.\/table_annovar.pl/\/docker_dependencies\/annovar\/table_annovar.pl/g' config.ini
RUN sed -i 's/.\/annotate_variation.pl/\/docker_dependencies\/annovar\/annotate_variation.pl/g' config.ini
ENV PATH "$PATH:/docker_dependencies/InterVar-2.2.1/"
WORKDIR "/"

RUN pip3 install pandas vcfpy biomart natsort pybedtools


WORKDIR "/"