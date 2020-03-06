# Make an image that has the dependencies for building GROMACS documentation.

# Make an intermediate image that can build a static Doxygen 1.8.5 that other
# containers will be able to use.

FROM ubuntu:18.04 as doxygen-builder
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /tmp
RUN \
  apt-get update && \
  apt-get -y -q=2 --no-install-suggests --no-install-recommends install \
      bison \
      build-essential \
      gcc \
      m4 \
      wget \
  && \
  wget --no-check-certificate https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/flex/2.5.35-10ubuntu3/flex_2.5.35.orig.tar.gz && \
  tar xf flex_2.5.35.orig.tar.gz && \
  cd flex-2.5.35 && \
  ./configure --prefix=/tmp/install-of-flex --disable-shared && \
  make -j && make install && cd .. && rm -rf flex* && \
  wget --no-check-certificate https://launchpad.net/ubuntu/+archive/primary/+sourcefiles/doxygen/1.8.5-1/doxygen_1.8.5.orig.tar.gz && \
  tar xf doxygen_1.8.5.orig.tar.gz && \
  cd doxygen-1.8.5 && \
  ./configure --flex /tmp/install-of-flex/bin/flex --static && \
  make -j && make install && cd .. && rm -rf doxygen* && \
  rm -rf /var/lib/apt/lists/*

# The ImageMagick package from apt has highly secure settings by
# default, suitable for use behind a webserver, which we don't
# need. So we use sed to remove those.
# We also install it separatly because it pulls in some dependencies
# that are needed for the documentation build.

FROM gromacs/ci-gcc-7:2020
WORKDIR /tmp
COPY --from=doxygen-builder /usr/local/bin/* /usr/local/bin/
RUN \
  apt-get update && \
  apt-get -y -q=2 --no-install-suggests --no-install-recommends install \
    graphviz \
    linkchecker \
    mscgen \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-fonts-recommended \
    texlive-fonts-extra && \
  apt-get -y install imagemagick && \
  rm -rf /var/lib/apt/lists/*
RUN \
  sed -i \
    '/\"XPS\"/d;/\"PDF\"/d;/\"PS\"/d;/\"EPS\"/d;/disable ghostscript format types/d' \
    /etc/ImageMagick-6/policy.xml && \
  pip3 install sphinx==1.6.1
