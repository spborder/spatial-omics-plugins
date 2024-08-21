# This Dockerfile is used to generate the docker image dsarchive/histomicstk
# This docker image includes the HistomicsTK python package along with its
# dependencies.
#
# All plugins of HistomicsTK should derive from this docker image

FROM rocker/r-ubuntu:20.04

LABEL maintainer="Sam Border CMI Lab <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get install --yes --no-install-recommends software-properties-common && \
    # As of 2018-04-16 this repo has the latest release of Python 2.7 (2.7.14) \
    # add-apt-repository ppa:jonathonf/python-2.7 && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get --yes --no-install-recommends -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libexpat1-dev \
    unzip \
    libhdf5-dev \
    libpython3-dev \
    python3.8-dev \
    python3.8-distutils \
    software-properties-common \
    libssl-dev \
    libffi-dev \
    # Standard build tools \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    # useful later \
    libmemcached-dev && \
    #apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN apt-get update ##[edited]

WORKDIR /
# Make Python3 the default and install pip.  Whichever is done last determines
# the default python version for pip.

#Make a specific version of python the default and install pip
RUN rm -f /usr/bin/python && \
    rm -f /usr/bin/python3 && \
    ln `which python3.8` /usr/bin/python && \
    ln `which python3.8` /usr/bin/python3 && \
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py && \
    rm get-pip.py
    #ln `which pip3` /usr/bin/pip 

RUN which  python && \
    python --version

RUN R -e 'install.packages("Seurat",repos ="https://cran.r-project.org")'
RUN R -e 'library(Seurat)'
ENV build_path=$PWD/build
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Copying codes files
ENV spat_path=$PWD
RUN mkdir -p $spat_path

RUN apt-get update && \
    apt-get install -y --no-install-recommends memcached && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY . $spat_path/
WORKDIR $spat_path

RUN pip install --no-cache-dir --upgrade --ignore-installed pip setuptools && \
    pip install --no-cache-dir .  && \
    rm -rf /root/.cache/pip/*

# Show what was installed
RUN python --version && pip --version && pip freeze

# define entrypoint through which all CLIs can be run
WORKDIR $spat_path/cli
LABEL entry_path = $spat_path/cli

# Test our entrypoint.  If we have incompatible versions of numpy and
# openslide, one of these will fail
RUN python -m slicer_cli_web.cli_list_entrypoint --list_cli
RUN python -m slicer_cli_web.cli_list_entrypoint SpotAnnotation --help

ENTRYPOINT ["/bin/bash", "docker-entrypoint.sh"]