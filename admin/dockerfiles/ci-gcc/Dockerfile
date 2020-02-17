# Make an image that has the dependencies for building GROMACS with gcc.


FROM gromacs/base:master
WORKDIR /tmp
ARG TOOL_VERSION
RUN \
  apt-get update && \
  apt-get -qqy --no-install-suggests --no-install-recommends install \
      gcc-$TOOL_VERSION \
      g++-$TOOL_VERSION && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/apt/archives/*

ENV CC=/usr/bin/gcc-$TOOL_VERSION
ENV CXX=/usr/bin/g++-$TOOL_VERSION

RUN \
  wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.2.tar.bz2 && \
  md5sum openmpi-4.0.2.tar.bz2 | grep 2fdfb4d31cfe1363d8bb5e7488e57cdf && \
  tar xvf openmpi-4.0.2.tar.bz2 && \
  (cd openmpi-4.0.2 && \
    ./configure \
        CC=$CC \
        CXX=$CXX && \
    make -j8 install) && \
  rm -rf openmpi-4.0.2*

# Alternatively, configure with --prefix=/usr
ENV LD_LIBRARY_PATH=/usr/local/lib
ENV MPICC=/usr/local/bin/mpicc
ENV MPICXX=/usr/local/bin/mpicxx
