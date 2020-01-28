# Provide an easy-to-reproduce environment in which to test full Python functionality.
# Produce an image with GROMACS installed. Use the root of the repository as the build context

# Optionally, set `--build-arg DOCKER_CORES=N` for a Docker engine running with access to more than 1 CPU.
#    REF=`git show -s --pretty=format:"%h"`
#    docker build -t gmxapi/gromacs:${REF} --build-arg DOCKER_CORES=4 -f gromacs.dockerfile ../..

# This image serves as a base for integration with the gmxapi Python tools and sample code.

ARG MPIFLAVOR=mpich
ARG REF=latest
FROM gmxapi/gromacs-dependencies-$MPIFLAVOR:$REF

ENV SRC_DIR /tmp/gromacs-source
COPY . $SRC_DIR

ENV BUILD_DIR /tmp/gromacs-build
RUN mkdir -p $BUILD_DIR
WORKDIR $BUILD_DIR

ARG DOCKER_CORES=1
# Allow the build type to be specified with `docker build --build-arg TYPE=something`
ARG TYPE=Release
RUN cmake $SRC_DIR \
        -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs \
        -DGMXAPI=ON \
        -DGMX_THREAD_MPI=ON \
        -DGMX_BUILD_HELP=OFF \
        -DGMX_REQUIRE_VALID_TOOLCHAIN=TRUE \
        -DCMAKE_BUILD_TYPE=$TYPE
RUN make -j$DOCKER_CORES
RUN make -j$DOCKER_CORES tests
RUN make -j$DOCKER_CORES install

# Default command provided for convenience since it inherits WORKDIR from above.
CMD make check
