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
    zip

RUN sed -i '/es_ES.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG es_ES.UTF-8
ENV LANGUAGE es_ES.UTF-8
ENV LC_ALL es_ES.UTF-8

RUN mkdir -p /docker_files
COPY ./docker_files/config_docker.json /docker_files

RUN mkdir -p /docker_dependencies
WORKDIR "/docker_dependencies"
ADD https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2 bcftools-1.20.tar.bz2
RUN tar -xf bcftools-1.20.tar.bz2
WORKDIR "/docker_dependencies/bcftools-1.20"
RUN ./configure --prefix=/
RUN make
RUN make install
RUN rm /docker_dependencies/bcftools-1.20.tar.bz2

WORKDIR "/docker_dependencies"
ADD https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz bedtools-2.31.1.tar.gz
RUN tar -xf bedtools-2.31.1.tar.gz
WORKDIR "/docker_dependencies/bedtools2"
RUN make
ENV PATH "$PATH:/docker_dependencies/bedtools2/bin/"
WORKDIR "/docker_dependencies"
RUN rm bedtools-2.31.1.tar.gz

COPY docker_files/annovar.latest.tar.gz /docker_dependencies/
WORKDIR "/docker_dependencies"
RUN tar xvzf annovar.latest.tar.gz
RUN rm annovar.latest.tar.gz


RUN mkdir -p /docker_directories/ref_genomes/37
RUN mkdir -p /docker_directories/ref_genomes/38
ADD https://zenodo.org/records/8045374/files/hs37d5.genome.tgz?download=1 /docker_directories/ref_genomes/37/hs37d5.genome.tgz
WORKDIR "/docker_directories/ref_genomes/37/"
RUN tar xvzf hs37d5.genome.tgz
RUN rm hs37d5.genome.tgz


WORKDIR "/docker_dependencies"
ADD https://github.com/WGLab/InterVar/archive/refs/tags/v2.2.1.tar.gz InterVar-2.2.1.tar.gz
RUN tar -xf InterVar-2.2.1.tar.gz
WORKDIR "/docker_dependencies/InterVar-2.2.1"
RUN sed -i 's/.\/convert2annovar.pl/\/docker_dependencies\/annovar\/convert2annovar.pl/g' config.ini
RUN sed -i 's/.\/table_annovar.pl/\/docker_dependencies\/annovar\/table_annovar.pl/g' config.ini
RUN sed -i 's/.\/annotate_variation.pl/\/docker_dependencies\/annovar\/annotate_variation.pl/g' config.ini
WORKDIR "/docker_dependencies/InterVar-2.2.1/"
RUN chmod -R 777 .
RUN chmod 755 Intervar.py
ENV PATH "$PATH:/docker_dependencies/InterVar-2.2.1/"
WORKDIR "/docker_dependencies"
RUN rm InterVar-2.2.1.tar.gz
WORKDIR "/"


WORKDIR "/docker_dependencies/InterVar-2.2.1/intervardb"
ADD https://zenodo.org/records/11177705/files/mim2gene.txt?download=1 mim2gene.txt
RUN chmod 755 mim2gene.txt
WORKDIR "/"
RUN mkdir -p "/docker_dependencies/InterVar-2.2.1/humandb"
WORKDIR "/docker_dependencies/InterVar-2.2.1/humandb"
RUN chmod -R 777 .
# Download databases to humandb
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene /docker_dependencies/InterVar-2.2.1//humandb/ 
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42a /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210501 /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_genome /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbscsnv11 /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb rmsk /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene /docker_dependencies/InterVar-2.2.1//humandb/
RUN perl /docker_dependencies/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene /docker_dependencies/InterVar-2.2.1//humandb/
WORKDIR "/"

RUN pip3 install pandas vcfpy biomart natsort pybedtools openpyxl --break-system-packages

RUN mkdir -p /docker_directories/categories/
COPY categories /docker_directories/categories

ADD https://zenodo.org/records/11146836/files/clinvar_database_GRCh37_20240421.txt?download=1 /docker_dependencies/clinvar/clinvar_database_GRCh37_20240421.txt
WORKDIR "/docker_dependencies/clinvar/"
RUN chmod 755 clinvar_database_GRCh37_20240421.txt

RUN mkdir -p /release_build/

WORKDIR "."
COPY modules /release_build/modules
COPY SFtool.py /release_build/SFtool.py
COPY docker_files/SFtool_singularity /release_build/SFtool
ENV PATH "$PATH:/release_build/"

WORKDIR "/"