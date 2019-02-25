# Update an ubuntu image with dependencies needed to build GROMACS and dependent packages.
# This version of the Dockerfile installs mpich.

# docker build -t gmxapi/gromacs-dependencies-mpich -f gromacs-dependencies.dockerfile .

# This image serves as a base for integration with the gmxapi Python tools and sample code.

FROM ubuntu:xenial

# Basic packages
RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install software-properties-common && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        git \
        libblas-dev \
        libcr-dev \
        libfftw3-dev \
        liblapack-dev \
        libxml2-dev \
        make \
        wget \
        zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# The CMake available with apt-get may be too old. Get v3.13.3 binary distribution.
ENV CMAKE_ROOT /usr/local/cmake
ENV PATH $CMAKE_ROOT/bin:$PATH

RUN set -ev && \
    mkdir /tmp/cmake && \
    (cd /tmp/cmake && \
    wget -c --no-check-certificate \
     https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-Linux-x86_64.tar.gz && \
    echo "78227de38d574d4d19093399fd4b40a4fb0a76cbfc4249783a969652ce515270  cmake-3.13.3-Linux-x86_64.tar.gz" \
     > cmake_sha.txt && \
    sha256sum -c cmake_sha.txt && \
    tar -xvf cmake-3.13.3-Linux-x86_64.tar.gz > /dev/null && \
    mv cmake-3.13.3-Linux-x86_64 $CMAKE_ROOT) && \
    rm -rf /tmp/cmake

RUN cmake --version

# mpich installation layer
RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        libmpich-dev \
        mpich && \
    rm -rf /var/lib/apt/lists/*
