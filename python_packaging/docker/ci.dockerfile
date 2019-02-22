# Provide an easy-to-reproduce environment in which to test full Python functionality.

# Run with default command or specify one of the scripts to be executed through the entrypoint.
#
#    docker run --rm -ti gmxapi/ci:0efcaed7b integrationtest.sh

# Building:
#
# Requires Docker 17.05 or higher.
#
# Optionally, set `--build-arg DOCKER_CORES=N` for a Docker engine running with access to more than 1 CPU.
#    REF=`git show -s --pretty=format:"%h"`
#    docker build -t gmxapi/ci:${REF} --build-arg REF=${REF} -f ci.dockerfile ..
#

#
# Use gromacs installation from gmxapi/gromacs image
#

ARG REF=master
FROM gmxapi/gromacs-mpich:$REF as gromacs
# This intermediate is necessary because the COPY command does not support syntax like the following:
#COPY --from=gmxapi/gromacs:$REF /usr/local/gromacs /usr/local/gromacs

FROM gmxapi/gromacs-dependencies-mpich

RUN apt-get update && \
    apt-get -yq --no-install-suggests --no-install-recommends install \
        python3 \
        python3-dev \
        python3-venv && \
    rm -rf /var/lib/apt/lists/*

# TODO: use pyenv for multiple Python installations.

ENV CMAKE_ROOT /usr/local/cmake
ENV PATH $CMAKE_ROOT/bin:$PATH

RUN groupadd -r testing && useradd -m -s /bin/bash -g testing testing

ADD --chown=testing:testing gmxapi/requirements-test.txt /home/testing/requirements.txt

USER testing

# TODO: Clean up pip cache.
RUN python3 -m venv $HOME/testing
RUN . $HOME/testing/bin/activate && \
    pip install --upgrade pip setuptools && \
    pip install jupyter
RUN . $HOME/testing/bin/activate && \
    pip install -r /home/testing/requirements.txt

COPY --from=gromacs $CMAKE_ROOT $CMAKE_ROOT
COPY --from=gromacs /usr/local/gromacs /usr/local/gromacs

ADD --chown=testing:testing gmxapi /home/testing/gmxapi
ADD --chown=testing:testing gmxapi/test /home/testing/gmxapi/test
ADD --chown=testing:testing gmxapi/src/gmxapi /home/testing/gmxapi/src/gmxapi

RUN . $HOME/testing/bin/activate && \
    . /usr/local/gromacs/bin/GMXRC && \
    (cd $HOME/gmxapi && \
     pip install . \
    )

ADD --chown=testing:testing acceptance /home/testing/acceptance
ADD --chown=testing:testing scripts /home/testing/scripts
ADD --chown=testing:testing test /home/testing/test

# TODO: this can be in the root user section above once it is stable
COPY docker/entrypoint.sh /

ENTRYPOINT ["/entrypoint.sh"]
CMD ["run_pytest"]


# MPI tests can be run in this container without requiring MPI on the host.
#     docker run --rm -t gmxapi/ci:${REF} /home/testing/scripts/run_pytest_mpi.sh
# We should also try tests with an MPI-connected set of docker containers.

# To be able to step through with gdb, run with something like the following, replacing
# 'imagename' with the name of the docker image built with this recipe.
# docker run --rm -ti --security-opt seccomp=unconfined imagename bash
