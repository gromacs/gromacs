# Make an image that has the basic dependencies for building GROMACS.
# This is the same for all other build images and gets used by those.

# Some optional GROMACS dependencies are obtained from the
# distribution, e.g.  fftw3, hwloc, blas and lapack so that the build
# is as fast as possible.
FROM ubuntu:18.04 as ci-basic-dependencies
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /tmp
RUN \
  apt-get update && \
  apt-get install -y \
    cmake \
    git \
    ninja-build \
    ccache \
    build-essential \
    wget \
    moreutils \
    rsync \
    libfftw3-dev \
    libhwloc-dev \
    liblapack-dev \
    xsltproc \
    python3-pip

# Make an image that has the dependencies for building GROMACS with gcc-7.
FROM ci-basic-dependencies as ci-gcc-7
WORKDIR /tmp
RUN apt-get -qqy --no-install-suggests --no-install-recommends install \
  gcc-7 \
  g++-7 && \
  rm -rf /var/lib/apt/lists/*
