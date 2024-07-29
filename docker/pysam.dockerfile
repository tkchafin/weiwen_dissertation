# Use an official Python runtime as a parent image
FROM python:3.12

# Set the working directory in the container
WORKDIR /app

# Install system dependencies and build tools
RUN apt-get update && apt-get install -y \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    bzip2 \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.17 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/samtools /usr/bin/samtools && \
    cd .. && rm -rf samtools-1.17 samtools.tar.bz2

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 -O bcftools.tar.bz2 && \
    tar -xjvf bcftools.tar.bz2 && \
    cd bcftools-1.17 && \
    make && \
    make prefix=/usr/local/bin install && \
    ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools && \
    cd .. && rm -rf bcftools-1.17 bcftools.tar.bz2

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.17 && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.17 htslib.tar.bz2

# Install Python dependencies
RUN pip install --no-cache-dir \
    pandas>2 \
    pysam

