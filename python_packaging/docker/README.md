# Docker images

This directory segregates the Dockerfiles to avoid clutter.

Assume you have already checked out the commit you want to build for.
Set an environment variable for convenience.

    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base gerrit_master HEAD`)
    REF=`git show -s --pretty=format:"%h"`

## Building

Note that the examples show me tagging the builds in the `gmxapi` dockerhub namespace.
If you don't have access to it, you can remove `gmxapi/` from the `-t` argument or use
a different dockerhub project space.

The different Dockerfiles require different build contexts.

For `gromacs-dependencies`, the build context doesn't matter.

    docker build -t gmxapi/gromacs-dependencies-mpich:${FORKPOINT} -f gromacs-dependencies.dockerfile .

For `gromacs`, the build context needs to be the root of the GROMACS repository.

    docker build -t gmxapi/gromacs-mpich:${REF} --build-arg DOCKER_CORES=4 --build-arg REF=${FORKPOINT} -f gromacs.dockerfile ../..

For integration testing here, we only want the `python_packaging` subdirectory.

    docker build -t gmxapi/ci-mpich:${REF} --build-arg REF=${REF} -f ci.dockerfile ..

## Running

### acceptance.dockerfile

For the interactive Jupyter notebook server, mapped to local port 8888:

    docker run --rm -t -p 8888:8888 gmxapi/acceptance

For a build that tests a feature branch that has been built and updated, add a tag, such as `gmxapi/acceptance:fr1`
    
For the a `bash` shell in a Conda environment under the `jovyan` user:

    docker run --rm -ti gmxapi/acceptance bash

### ci.dockerfile

The `entrypoint.sh` script activates the python venv and wraps commands in a `bash` `exec`.
The default command is a script sourced from `../scripts/run_pytest.sh`. You can use this,
other scripts, `bash`, etc.

    docker run --rm -t gmxapi/ci-mpich:${REF}
    docker run --rm -t gmxapi/ci-mpich:${REF} run_pytest_mpi
    docker run --rm -t gmxapi/ci-mpich:${REF} bash

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

Travis-CI builds and pushes a chain of Docker images to the `gmxapi` dockerhub organization.
The `kassonlab` GitHub organization `gromacs-gmxapi` repository branches that are descended from the `kassonLabFork`
branch have the necessary Travis-CI configuration.
