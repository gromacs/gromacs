# Provide an easy-to-reproduce environment in which to test full Python functionality.

# Run with default command or specify one of the scripts to be executed through the entrypoint.
#
#    docker run --rm -ti gmxapi/ci-mpich:fr1 integrationtest
#
# Building:
#
# Requires Docker 17.05 or higher.
#
# Note to maintainers:
# Build from the GROMACS image at the current fork point. Tag with the feature
# name or the current revision.
#
#    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base master HEAD`)
#    REF=`git show -s --pretty=format:"%h"`
#    # or
#    REF="fr1"
#    docker build -t gmxapi/ci-mpich:${REF} --build-arg REF=${FORKPOINT} -f ci.dockerfile ..
#

ARG REF=latest
ARG MPIFLAVOR=mpich

FROM gmxapi/gromacs-dependencies-$MPIFLAVOR as python-base

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -yq --no-install-suggests --no-install-recommends \
    install \
        python3 \
        python3-dev \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

# TODO: Use non-system Python installations for explicit version coverage.
# Consider pyenv for generic management of Python environment.

RUN groupadd -r testing && useradd -m -s /bin/bash -g testing testing

USER testing

ENV VENV /home/testing/venv
RUN python3 -m venv $VENV
RUN . $VENV/bin/activate && \
    pip install --no-cache-dir --upgrade pip setuptools wheel

ADD --chown=testing:testing requirements-*.txt /home/testing/gmxapi/

#
# Use gromacs installation from gmxapi/gromacs image
#

FROM gmxapi/gromacs-$MPIFLAVOR:$REF as gromacs
# This intermediate is necessary because the COPY command does not support syntax like the following:
#COPY --from=gmxapi/gromacs:$REF /usr/local/gromacs /usr/local/gromacs

FROM python-base

COPY --from=gromacs /usr/local/gromacs /usr/local/gromacs

RUN $VENV/bin/python -m pip install --upgrade pip setuptools wheel
RUN $VENV/bin/python -m pip install --no-cache-dir --no-build-isolation -r /home/testing/gmxapi/requirements-test.txt

ADD --chown=testing:testing src /home/testing/gmxapi/src
ADD --chown=testing:testing src/gmxapi /home/testing/gmxapi/src/gmxapi

# We use "--no-cache-dir" to reduce Docker image size.
RUN . $VENV/bin/activate && \
    (cd $HOME/gmxapi/src && \
     CMAKE_ARGS="-Dgmxapi_ROOT=/usr/local/gromacs -C /usr/local/gromacs/share/cmake/gromacs/gromacs-hints.cmake" \
      pip install --no-cache-dir --verbose . \
    )

ADD --chown=testing:testing sample_restraint /home/testing/sample_restraint

# To test behavior as in GitLab CI, copy the googletest sources, export CI=1 to the cmake
# configure command, and remove the option to download googletest.
#COPY --from=gromacs /gromacs-source/src/external/googletest /home/testing/sample_restraint/external/googletest
RUN cmake --version
RUN . $VENV/bin/activate && \
    . /usr/local/gromacs/bin/GMXRC && \
    (cd $HOME/sample_restraint && \
     mkdir build && \
     cd build && \
     cmake .. \
             -DPYTHON_EXECUTABLE=$VENV/bin/python \
             -DDOWNLOAD_GOOGLETEST=ON \
             -DGMXAPI_EXTENSION_DOWNLOAD_PYBIND=ON && \
     make -j4 && \
     make tests && \
     make test && \
     make install \
    )

ADD --chown=testing:testing src/test /home/testing/gmxapi/test
ADD scripts /docker_entry_points
COPY docker/entrypoint.sh /

ENTRYPOINT ["/entrypoint.sh"]
CMD ["run_full"]


# MPI tests can be run in this container without requiring MPI on the host.
# (We suggest running your docker engine with multiple CPU cores allocated.)
#     docker run --rm -t gmxapi/ci:${REF} /home/testing/scripts/run_full_mpi.sh
# We should also try tests with an MPI-connected set of docker containers.

# To be able to step through with gdb, run with something like the following, replacing
# 'imagename' with the name of the docker image built with this recipe.
# docker run --rm -ti --security-opt seccomp=unconfined imagename bash
