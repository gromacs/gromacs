# Make a base container useful for all others

FROM ubuntu:18.04 as base
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update

# Make a container that can build a static Doxygen 1.8.5 that other
# containers will be able to use.

FROM base as doxygen-builder
RUN apt-get install -y \
  gcc \
  build-essential \
  m4 \
  bison \
  wget
WORKDIR /tmp
RUN wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/flex/2.5.35-10ubuntu3/flex_2.5.35.orig.tar.gz && \
  tar xf flex_2.5.35.orig.tar.gz && \
  cd flex-2.5.35 && \
  ./configure --prefix=/tmp/install-of-flex --disable-shared && \
  make -j && make install && cd .. && rm -rf flex*
RUN wget https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/doxygen/1.8.5-1/doxygen_1.8.5.orig.tar.gz && \
  tar xf doxygen_1.8.5.orig.tar.gz && \
  cd doxygen-1.8.5 && \
  ./configure --flex /tmp/install-of-flex/bin/flex --static && \
  make -j && make install && cd .. && rm -rf doxygen*

# Make an image that has the basic dependencies for building GROMACS.

# TODO Remove git dependency?

FROM base as ci-basic-dependencies
WORKDIR /tmp
RUN apt-get install -y \
  cmake \
  git \
  ninja-build \
  ccache

# Make an image that has the dependencies for building GROMACS documentation.

# The ImageMagick package from apt has highly secure settings by
# default, suitable for use behind a webserver, which we don't
# need. So we use sed to remove those.

# Some optional GROMACS dependencies are obtained from the
# distribution, e.g.  fftw3, hwloc, blas and lapack so that the build
# is as fast as possible.

FROM ci-basic-dependencies as ci-docs
WORKDIR /tmp
COPY --from=doxygen-builder /usr/local/bin/* /usr/local/bin/
RUN apt-get install -y \
  texlive-latex-base \
  texlive-fonts-recommended \
  texlive-fonts-extra \
  texlive-latex-extra
RUN \
  apt-get install -y \
    clang \
    graphviz \
    imagemagick \
    linkchecker \
    mscgen \
    liblapack-dev \
    libfftw3-dev \
    libhwloc-dev \
    python3-pip && \
  sed -i \
    '/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d' \
    /etc/ImageMagick-6/policy.xml && \
  pip3 install sphinx==1.6.1
