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

# Make an image that has the dependencies for building GROMACS with clang-7.
FROM ci-basic-dependencies as ci-clang-7
WORKDIR /tmp
RUN \
  apt-get -y -q=2 --no-install-suggests --no-install-recommends install \
    apt-utils \
    software-properties-common \
    gpg-agent && \
  wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - && \
  apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-7 main" && \
  apt-get -qq update && \
  apt-get -qqy --no-install-suggests --no-install-recommends install \
    clang++-7 \
    clang-tools-7 \
    libomp-dev && \
  rm -rf /var/lib/apt/lists/*

# Make a container that can build a static Doxygen 1.8.5 that other
# containers will be able to use.

FROM ubuntu:18.04 as doxygen-builder
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /tmp
RUN apt-get install -y \
  gcc \
  build-essential \
  m4 \
  bison \
  wget && \
  wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/flex/2.5.35-10ubuntu3/flex_2.5.35.orig.tar.gz && \
  tar xf flex_2.5.35.orig.tar.gz && \
  cd flex-2.5.35 && \
  ./configure --prefix=/tmp/install-of-flex --disable-shared && \
  make -j && make install && cd .. && rm -rf flex* && \
  wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/doxygen/1.8.5-1/doxygen_1.8.5.orig.tar.gz && \
  tar xf doxygen_1.8.5.orig.tar.gz && \
  cd doxygen-1.8.5 && \
  ./configure --flex /tmp/install-of-flex/bin/flex --static && \
  make -j && make install && cd .. && rm -rf doxygen*

# Make an image that has the dependencies for building GROMACS documentation.

# The ImageMagick package from apt has highly secure settings by
# default, suitable for use behind a webserver, which we don't
# need. So we use sed to remove those.

# Some optional GROMACS dependencies are obtained from the
# distribution, e.g.  fftw3, hwloc, blas and lapack so that the build
# is as fast as possible.

FROM ci-clang-7 as ci-docs
WORKDIR /tmp
COPY --from=doxygen-builder /usr/local/bin/* /usr/local/bin/
RUN apt-get install -y \
  texlive-latex-base \
  texlive-fonts-recommended \
  texlive-fonts-extra \
  texlive-latex-extra &&
  apt-get install -y \
    graphviz \
    imagemagick \
    linkchecker \
    mscgen && \
  sed -i \
    '/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d' \
    /etc/ImageMagick-6/policy.xml && \
  pip3 install sphinx==1.6.1 && \
  rm -rf /var/lib/apt/lists/*
