# Make an image that has the basic dependencies for building GROMACS.
# This is the same for all other build images and gets used by those.

# Some optional GROMACS dependencies are obtained from the
# distribution, e.g.  fftw3, hwloc, blas and lapack so that the build
# is as fast as possible.
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /tmp
RUN \
  apt-get update && \
  apt-get -y -q=2 --no-install-suggests --no-install-recommends install \
    build-essential \
    ccache \
    cmake \
    git \
    libfftw3-dev \
    libhwloc-dev \
    liblapack-dev \
    moreutils \
    ninja-build \
    python3-pip \
    rsync \
    wget \
    xsltproc \
    && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/apt/archives/*
