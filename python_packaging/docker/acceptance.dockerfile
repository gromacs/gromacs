# Provide an easy-to-reproduce environment in which to run the exploratory
# notebook or acceptance tests.
#
# Build GROMACS with modifications for gmxapi Python package testing.
#     docker build -t acceptance --build-arg DOCKER_CORES=4 -f acceptance.dockerfile ..
#
# Note that the build context is the `python_packaging` subdirectory, not the root
# of the repository. The root `.dockerignore` file excludes the `python_packagin`
# subdirectory to avoid overly conservative expiration of docker build caches.
#
# Note to maintainers: Update public feature-test images
#     docker tag acceptance gmxapi/acceptance:fr1
#     docker push gmxapi/acceptance:fr1
#
# Either run for an interactive notebook server or in various unattended modes.
# For jupyter notebook server mapped to local port 8888:
#
#    docker run --rm -ti -p 8888:8888 gmxapi/acceptance
#
# For bash shell as the "jovyan" user:
#
#    docker run --rm -ti gmxapi/acceptance bash

FROM jupyter/scipy-notebook

USER root
RUN \
    apt-get update && \
    apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        gdb \
        git \
        libmpich-dev \
        mpich \
        vim \
        ssh \
        wget && \
    rm -rf /var/lib/apt/lists/*

USER jovyan
RUN conda install pytest tox cmake

# Note: This is a little sloppy, and assumes gmxapi/gromacs-mpich:latest is
# based on a compatible base image.
# TODO: Split out a acceptance-dependencies.dockerfile? Accept slower builds and larger images?
COPY --from=gmxapi/gromacs-mpich /usr/local/gromacs /usr/local/gromacs

#
# Get RequiredFunctionality notebook
#

COPY acceptance /home/jovyan/acceptance

#
# gmxapi Python package layer
#

ADD --chown=jovyan:users gmxapi /home/jovyan/gmxapi
RUN . /usr/local/gromacs/bin/GMXRC && \
    cd /home/jovyan/gmxapi && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir .

# MPI tests can be run in this container without requiring MPI on the host.
# We should also try tests with an MPI-connected set of docker containers.

# To be able to step through with gdb, run with something like the following, replacing
# 'imagename' with the name of the docker image built with this recipe.
# docker run --rm -ti --security-opt seccomp=unconfined imagename bash
