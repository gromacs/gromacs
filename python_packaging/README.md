# Python package sources

This project attempts to update the gmxapi efforts from GROMACS 2019,
merge external repositories from
https://github.com/kassonlab/gmxapi
and
https://github.com/kassonlab/sample_restraint,
and build the new functionality proposed at
https://github.com/kassonlab/gmxapi-scripts

This directory provides the files that will be copied to the GROMACS installation location from which users may 
install Python packages.
This allows C++ extension modules to be built against a user-chosen GROMACS installation,
but for a Python interpreter that is very likely different from that used
by the system administrator who installed GROMACS.

## Repository organization

https://github.com/kassonlab/gromacs-gmxapi has a branch called `kassonLabFork`
that is intended as a patch to the GROMACS master branch to allow separable, testable
gmxapi feature development. The branch adds a `.travis.yml` and other differences
supporting development and testing by the Kasson Lab.

A squash of this branch could be used as the basis of a sandbox branch on Gerrit
and the branch could be re-forked from there.

GROMACS `master` can be merged into `kassonLabFork` as necessary or as feature
requests are resolved.

The `kassonlab` GitHub repository has feature-development branches from which Gerrit
changes are intended to be extracted. Gerrit patch sets should be equivalent to the
`diff` applied by the feature branch to the `kassonLabFork` branch. In other words,
the Gerrit patch set for a proposed resolution to feature request `fr1`
is derived from.

    git rebase --onto gerrit_master kassonLabFork fr1

## gmxapi

Python framework for gmxapi high-level interface. Contains functionality not ready or not destined for gromacs 
bindings package.

Current packaging:

    # Install GROMACS to /path/to/gromacs
    source /path/to/gromacs/bin/GMXRC
    source $HOME/somevirtualenv/bin/activate
    (cd src && pip install -r requirements.txt && pip install .)
    python -c 'import gmxapi as gmx'
    pytest src/test
    pytest test

Intended packaging:

    # Build and install GROMACS, which also builds and installs a Python sdist
    source $HOME/somevirtualenv/bin/activate
    pip install -r /path/to/gromacs/python-requirements.txt
    pip install /path/to/gromacs/gmxapi*tar.gz
    python -c 'import gmxapi as gmx'

We can also consider PyPI source and binary distributions in the future.

## Docker and Travis-CI testing

The Python packaging will be tested on Jenkins with Docker-based CI, but this
infrastructure is a little way off. In the mean time, we are trying to submit
changes that do not affect main line GROMACS development or building and to
perform testing with Docker externally. Users may build and run Docker images
from the `python_packaging/docker` directory. The `kassonlab` Travis-CI account
will automatically build and run Docker-based tests for each outstanding
feature branch.

The Dockerfiles direct a few different Linux, Python, and GROMACS configurations,
build and install the `gmxapi` and `sample_restraint` packages, and provide a
few styles of testing through the `scripts` accessible through `entrypoint.sh`.

In successive build stages, Travis-CI is directed to use a series of Docker images,
referred to here with their dockerhub repository, an explanation of tags,
and the Dockerfiles from which they are built.

1. `gmxapi/gromacs-dependencies-<matrix>:<tag>` Ubuntu distribution with dependencies for
   various GROMACS build configurations. `<matrix>` encodes the build matrix dimensions
   for things like compiler and MPI flavor. Travis-CI will rebuild this for commits to
   `kassonLabFork`, but if the `<matrix>` string changes, a more privileged dockerhub
   account will have to push the new repository the first time and then grant access
   to the service account used by the Travis-CI configuration.
   `<tag>` is the (short) git revision hash
   of the `master` branch commit corresponding to the current state of the `kassonLabFork`
   branch.
2. `gmxapi/gromacs-<matrix>:<tag>` Builds on `gromacs-dependencies-<matrix>`, where
   `<matrix>` has the same meaning as above. `<tag>` is the (short) git revision hash
   of the `master` branch commit corresponding to the current state of the `kassonLabFork`
   branch.
   This is recorded in the `kassonLabFork` `.travis.yml`.
3. `gmxapi/ci-<matrix>:<tag>` starts with `gromacs-<matrix>` and merges in the
    `python_packaging` changes associated with the feature branch indicated by `<tag>`

`acceptance.dockerfile` is based on the `jupyter/scipy-notebook` image and is probably
too expensive to build in Travis-CI builds right now,
but will be available to pull from `gmxapi/acceptance:frN` as
feature test builds are available.
It is primarily intended for running the acceptance tests interactively,
though it also provides a useful alternative environment,
in that the base Python installation is a user-space Conda install.

Additional information in `python_packaging/docker/README.md`.

Hint: the fork point from `master` and the current git ref can be set as environment variables:

    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base gerrit_master HEAD`)
    REF=`git show -s --pretty=format:"%h"`

## pybind11

Python bindings are expressed in C++ using the pybind11 template library.
The pybind11 repository is mirrored in GROMACS project sources and
installed with GROMACS for convenience and reproducibility.

pybind11 sources were added to the repository with `git-subtree` and a squash merge.

To update the version of pybind11 included with the GROMACS repository,
`git-subtree` uses meta-data stashed in the git history to prepare an appropriate squash merge.

# Cross compiling

On some systems, GROMACS will have been built and installed for a different
architecture than the system on which the Python package will be compiled.
We need to use CMake Tool-chains to support cross-compiling for the target architecture.

Note: scikit-build can use CMake Toolchains to properly handle `pip` builds.

# Offline installation

The `pip install` options `--no-index` and `--find-links` allow for an offline stash of package archives so that 
satisfying dependencies for a new virtualenv does not require network access or lengthy build times.
