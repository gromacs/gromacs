# Make an image that has the dependencies for building GROMACS with gcc.


FROM gromacs/base:2020
WORKDIR /tmp
ARG TOOL_VERSION
RUN \
  apt-get update && \
  apt-get -qqy --no-install-suggests --no-install-recommends install \
      gcc-$TOOL_VERSION \
      g++-$TOOL_VERSION && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/apt/archives/*
