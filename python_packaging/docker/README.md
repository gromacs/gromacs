# Docker images

This directory segregates the Dockerfiles to avoid clutter. The Dockerfiles
here help to build and test gmxapi software. They may be subsumed or supplanted
by future infrastructure.

Assume you have already checked out the commit you want to build for.
Assume the following definitions.

    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base master HEAD`)
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

**Warning:** The `--rm` flag tells Docker to remove the container after the
process completes. This prevents unnamed containers from consuming more disk
space on the host with each `run`.

Alternatively, replace `--rm` with `--name containername` to save a named copy
of the container and consider using `commit`ing snapshots of the container.

Refer to Docker documentation for details.

### ci.dockerfile

The default user for this image is `testing`. The gmxapi and sample_restraint
Python packages are installed into a Python 3 `venv` (virtual environment) owned
by the `testing` user.

The `entrypoint.sh` script activates the python venv and wraps commands in a `bash` `exec`.
The default command is a script sourced from `../scripts/run_pytest.sh`. You can use this,
other scripts, `bash`, etc.

    docker run --rm -t gmxapi/ci-mpich:${TAG}
    docker run --rm -t gmxapi/ci-mpich:${TAG} run_pytest_mpi
    docker run --rm -ti gmxapi/ci-mpich:${TAG} bash

### Why venv?

`venv` is the suggested and primary installation mode for the gmxapi Python package,
so it is the most important installation mode to test.

These scripts will ultimately be ported to as-yet-undefined GROMACS testing
infrastructure.
It seems equally plausible to have a single image with multiple Python installations
as to have multiple Docker images, or that single-Python docker images would use
some sort of Python virtualenv system to manage non-default Python interpreters.

Since the installation, tests, and shell environment for post-mortem all use the
"testing" user instead of "root", the venv provides a tidy place to work, avoids
the need for the `--user` flag to `pip`, and gives us a convenient place to do
`pip freeze` to get a known baseline of a working Python environment.

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

### notebook.dockerfile

Built on `gmxapi/ci-mpich:latest`, this image adds a `notebook` entry point to
be used as the new default command. Run with port mapping for http port 8888 for
easy access to a Jupyter notebook server running in the docker container's
`testing` user home directory.

    docker run --rm -ti -p 8888:8888 gmxapi/notebook

Note that, when run with `--rm`, changes made to files in the container will be
lost when `docker run` completes.
To preserve work in a notebook run this way,
download the `ipynb` through theJupyter web interface
(such as when updating the examples in the repository).

### docs.dockerfile

For very quick and isolated documentation builds on top of the gmxapi/ci-mpich 
image, build the image from docs .dockerfile.
The resulting image is a small web server image (without GROMACS or gmxapi installed) 
with html content built in and copied from a temporary container.

    docker run --rm -p 8080:80 gmxapi/docs

Then browse to http://localhost:8080/

## Automation

*TODO: Update this section as CI infrastructure evolves.*

Travis-CI builds and pushes a chain of Docker images to the `gmxapi` dockerhub organization.
The `kassonlab` GitHub organization `gromacs-gmxapi` repository branches that are descended from the `kassonLabFork`
branch have the necessary Travis-CI configuration.
