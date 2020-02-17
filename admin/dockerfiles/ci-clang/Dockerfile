# Make an image that has the dependencies for building GROMACS with clang.
# Note when specifying TOOL_VERSION that clang 6.0 packages use the minor version
# in the name, while 7 and 8 do not.
FROM gromacs/base:master
WORKDIR /tmp
ARG TOOL_VERSION
RUN \
  apt-get update && \
  apt_version=$TOOL_VERSION && \
  if [ "$TOOL_VERSION" -lt "7" ] ; then apt_version="$apt_version.0"; fi && \
  apt-get -y -q=2 --no-install-suggests --no-install-recommends install \
    apt-utils \
    software-properties-common \
    gpg-agent && \
  wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - && \
  apt-add-repository "deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-$apt_version main" && \
  apt-get -qq update && \
  apt-get -qqy --no-install-suggests --no-install-recommends install \
    clang++-$apt_version \
    clang-tools-$apt_version \
    libomp-dev && \
  if [ "$apt_version" != "$TOOL_VERSION" ] ; then \
    ln -s /usr/bin/clang-$apt_version /usr/bin/clang-$TOOL_VERSION; \
    ln -s /usr/bin/clang++-$apt_version /usr/bin/clang++-$TOOL_VERSION; \
    ln -s /usr/bin/clang-format-$apt_version /usr/bin/clang-format-$TOOL_VERSION; \
    ln -s /usr/bin/clang-tidy-$apt_version /usr/bin/clang-tidy-$TOOL_VERSION; fi && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/apt/archives/*

ENV CC=/usr/bin/clang-$TOOL_VERSION
ENV CXX=/usr/bin/clang++-$TOOL_VERSION

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
