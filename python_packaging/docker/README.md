# Docker images

This directory segregates the Dockerfiles to avoid clutter. The Dockerfiles
here help to build and test gmxapi software. They may be subsumed or supplanted
by future Jenkins infrastructure.

Assume you have already checked out the commit you want to build for.
Assume the following definitions.

    git fetch https://github.com/gromacs/gromacs.git master
    git branch gerrit_master FETCH_HEAD
    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base gerrit_master HEAD`)
    TAG="fr1" # for functional requirement 1

## Building

Note that the examples show the builds tagged in the `gmxapi` dockerhub namespace.
If you don't have access to it, you can remove `gmxapi/` from the `-t` argument or use
a different dockerhub project space.

The different Dockerfiles require different build contexts (final path argument).

For `gromacs-dependencies`, the build context doesn't matter. Just use `.` in the
`docker` directory.

    docker build -t gmxapi/gromacs-dependencies-mpich -f gromacs-dependencies.dockerfile .
    # optional:
    docker tag gmxapi/gromacs-dependencies-mpich gmxapi/gromacs-dependencies-mpich:${FORKPOINT}

This should rarely be necessary, and the dependent images can probably just pull the `latest`
from dockerhub.

For `gromacs`, the build context needs to be the root of the GROMACS repository (`../..`).
In case images for feature branches diverge too much or become tightly coupled to particular revisions in `master`,
it may be useful to tag this image to annotate the GROMACS build.

    # Use DOCKER_CORES to let `make` use all cores available to the Docker engine.
    # optionally include an additional `--build-arg REF=${FORKPOINT}`
    docker build -t gmxapi/gromacs-mpich --build-arg DOCKER_CORES=4 -f gromacs.dockerfile ../..
    # optional:
    docker tag gmxapi/gromacs-mpich gmxapi/gromacs-mpich:${FORKPOINT}

For integration testing here, we only want the `python_packaging` subdirectory (`..`).
The image should be tagged according to the functionality it is intended to demonstrate.

    docker build -t gmxapi/ci-mpich:${TAG} --build-arg REF=${FORKPOINT} -f ci.dockerfile ..

## Running

### ci.dockerfile

The `entrypoint.sh` script activates the python venv and wraps commands in a `bash` `exec`.
The default command is a script sourced from `../scripts/run_pytest.sh`. You can use this,
other scripts, `bash`, etc.

    docker run --rm -t gmxapi/ci-mpich:${TAG}
    docker run --rm -t gmxapi/ci-mpich:${TAG} run_pytest_mpi
    docker run --rm -ti gmxapi/ci-mpich:${TAG} bash

#### Entry points

Images have an ENTRYPOINT script that allows
a container to be launched with a reasonable default command, an executable
locatable in the container, or one of the scripts copied from the `scripts`
directory parallel to this one. These additional scripts are primarily to
allow easy execution of certain suites of tests.

See the `scripts` directory for details.

### Debugging

To be able to step through with gdb, run with something like the following, replacing
'imagename' with the name of the docker image built with this recipe.

    docker run --rm -ti --security-opt seccomp=unconfined imagename bash

## Automation

*TODO: Update this section as Jenkins infrastructure evolves.*

Travis-CI builds and pushes a chain of Docker images to the `gmxapi` dockerhub organization.
The `kassonlab` GitHub organization `gromacs-gmxapi` repository branches that are descended from the `kassonLabFork`
branch have the necessary Travis-CI configuration.
