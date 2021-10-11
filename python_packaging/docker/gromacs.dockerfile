# Provide an easy-to-reproduce environment in which to test full Python functionality.
# Produce an image with GROMACS installed. Use the root of the repository as the build context

# Optionally, set `--build-arg DOCKER_CORES=N` for a Docker engine running with access to more than 1 CPU.
#    REF=`git show -s --pretty=format:"%h"`
#    docker build -t gmxapi/gromacs-${MPIFLAVOR}:${REF} \
#               --build-arg DOCKER_CORES=4 \
#               --build-arg MPIFLAVOR=${MPIFLAVOR} \
#               -f gromacs.dockerfile ../..

# This image serves as a base for integration with the gmxapi Python tools and sample code.

ARG MPIFLAVOR=mpich
ARG REF=latest
FROM gmxapi/gromacs-dependencies-$MPIFLAVOR:$REF as build

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
        -DGMX_USE_RDTSCP=OFF \
        -DGMX_INSTALL_LEGACY_API=ON \
        -DCMAKE_BUILD_TYPE=$TYPE
RUN make -j$DOCKER_CORES
RUN make -j$DOCKER_CORES tests
RUN make -j$DOCKER_CORES install

FROM gmxapi/gromacs-dependencies-$MPIFLAVOR:$REF

COPY --from=build /usr/local/gromacs /usr/local/gromacs
