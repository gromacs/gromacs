# Update an ubuntu image with dependencies needed to build GROMACS and dependent packages.
# This version of the Dockerfile installs mpich.

# docker build -t gmxapi/gromacs-dependencies-mpich -f gromacs-dependencies.dockerfile .

# This image serves as a base for integration with the gmxapi Python tools and sample code.

FROM ubuntu:focal as base

# Basic packages
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        software-properties-common \
        build-essential && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        cmake \
        git \
        libblas-dev \
        libfftw3-dev \
        liblapack-dev \
        libxml2-dev \
        make \
        vim \
        wget \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# mpich installation layer
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        libmpich-dev \
        mpich && \
    rm -rf /var/lib/apt/lists/*
