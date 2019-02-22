# Provide an easy-to-reproduce environment in which to run the exploratory
# notebook or acceptance tests.
#
# Build GROMACS with modifications for gmxapi Python package testing.
#
# Either run for an interactive notebook server or in various unattended modes.
# For jupyter notebook server mapped to local port 8888:
#
#    docker run --rm -ti -p 8888:8888 gmxapi/acceptance
#
# For bash shell as the "jovyan" user:
#
#    docker run --rm -ti gmxapi/acceptance bash
#
# TODO: instructions and automation for running a Python script extracted from the IPython notebook.

FROM jupyter/scipy-notebook as intermediate

# jupyter/scipy-notebook sets the USER to jovyan. Temporarily switch it back to root to install more stuff.
USER root

#
# Dependencies layer
#

RUN \
    apt-get update && \
    apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        byobu \
        curl \
        libfftw3-dev \
        libopenblas-dev \
        libopenmpi-dev \
        gdb \
        git \
        graphviz \
        htop \
        man \
        mscgen \
        openmpi-bin \
        ssh \
        software-properties-common \
        unzip \
        vim \
        wget && \
    rm -rf /var/lib/apt/lists/*
# By putting the above commands together, we avoid creating a large intermediate layer with all of /var/lib/apt/lists

#
# CMake layer
#
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

#
# GROMACS layer
#

ENV SRC_DIR /tmp/gromacs-source
COPY . $SRC_DIR

ENV BUILD_DIR /tmp/gromacs-build
RUN mkdir -p $BUILD_DIR

ARG DOCKER_CORES=1
# Allow the build type to be specified with `docker build --build-arg TYPE=something`
ARG TYPE=Release
# Note: AVX2 instructions not available in older docker engines.
RUN cd $BUILD_DIR && \
    cmake $SRC_DIR \
        -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs \
        -DGMXAPI=ON \
        -DGMX_THREAD_MPI=ON \
        -DGMX_BUILD_HELP=OFF \
        -DGMX_SIMD=AVX_256 \
        -DGMX_USE_RDTSCP=OFF \
        -DGMX_HWLOC=OFF \
        -DCMAKE_BUILD_TYPE=$TYPE && \
    make -j$DOCKER_CORES install

USER jovyan
RUN conda install pytest tox cmake

ADD --chown=jovyan:users python_packaging/acceptance /home/jovyan/acceptance

# MPI tests can be run in this container without requiring MPI on the host.
# We should also try tests with an MPI-connected set of docker containers.

# To be able to step through with gdb, run with something like the following, replacing
# 'imagename' with the name of the docker image built with this recipe.
# docker run --rm -ti --security-opt seccomp=unconfined imagename bash
