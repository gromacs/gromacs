# Python package sources

This directory exists as a staging area supporting GROMACS enhancement
[#2045](https://redmine.gromacs.org/issues/2045),
which attempts to update the gmxapi efforts from GROMACS 2019,
merge external repositories from
https://github.com/kassonlab/gmxapi
and
https://github.com/kassonlab/sample_restraint,
and build the new functionality proposed at
https://github.com/kassonlab/gmxapi-scripts

## Repository organization

**TODO: testing infrastructure and project management conventions to allow fully integrated development**
**TODO: Consider long term homes of these directory contents.**

## gmxapi

Python framework for gmxapi high-level interface.

The `src` directory provides the files that will be copied to the GROMACS installation location from which users may 
install Python packages.
This allows C++ extension modules to be built against a user-chosen GROMACS installation,
but for a Python interpreter that is very likely different from that used
by the system administrator who installed GROMACS.

To build and install the Python package,
first install GROMACS to `/path/to/gromacs`.
Then, install the package in a Python virtualenv.

    source /path/to/gromacs/bin/GMXRC
    python3 -m venv $HOME/somevirtualenv
    source $HOME/somevirtualenv/bin/activate
    (cd src && pip install -r requirements.txt && pip install .)
    python -c 'import gmxapi as gmx'

Use `pytest` to run unit tests and integration tests.

    pip install -r requirements-test.txt
    pytest src/test
    pytest test

For additional discussion on packaging and distribution, see
https://redmine.gromacs.org/issues/2896

## Sample MD extension code

`sample_restraint` is a subtree containing a complete CMake project for building
pluggable GROMACS MD extensions for execution through gmxapi. Up to and
including version 0.0.7 of the sample code, the sub-project lived at
https://github.com/kassonlab/sample_restraint/ and was supported by GROMACS 2019.

The GROMACS repository becomes the upstream source for the sample code for
GROMACS releases 2020 and higher. Refer to [README.md](sample_restraint/README.md)
in `python_packaging/sample_restraint` for more information.

To use a plugin based on the sample restraint, you will need to build and install
the gmxapi Python package (above).

**todo** CookieCutter repo for easy forking.

## Docker and Travis-CI testing

**TODO: Migrate to Jenkins-based CI as Docker infrastructure becomes available.**
Infastructure described here is transitional and reflects our need to be able to see code work 
in order to review it satisfactorily in the period before GROMACS CI infrastructure 
is ready for the load. At some point the Docker aspects will change,
or be removed as appropriate.

The Python packaging will be tested on Jenkins with Docker-based CI, but this
infrastructure is a little way off. In the mean time, we are trying to submit
changes that do not affect main line GROMACS development or building and to
perform testing with Docker externally. Users may build and run Docker images
from the `python_packaging/docker` directory. The `kassonlab` Travis-CI account
will automatically build and run Docker-based tests for each outstanding
feature branch.

The Dockerfiles
* direct a few different Linux, Python, and GROMACS configurations,
* build and install the `gmxapi` and `sample_restraint` packages, and
* provide a few styles of testing through the `scripts` accessible through `entrypoint.sh`.

In successive build stages, Travis-CI is directed to use a series of Docker images,
referred to here with their dockerhub repository, an explanation of tags,
and the Dockerfiles from which they are built.
The image naming scheme encodes a build matrix element in the repository name and
a git branch or reference in the image tag (described below).
Additional information in `python_packaging/docker/README.md`.

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

Hint: the fork point from `master` and the current git ref can be set as environment variables:

    FORKPOINT=$(git show -s --pretty=format:"%h" `git merge-base gerrit_master HEAD`)
    REF=`git show -s --pretty=format:"%h"`

## External project code

Refer to `./src/external/README.md` for current details on the copied external
sources.

# scikit-build

For the C++ extension module `_gmxapi`, 
[scikit-build](https://scikit-build.readthedocs.io/en/latest/)
provides glue between Python package tools and CMake infrastructure in `./src`.
Scikit-build is installed with Python packaging tools automatically with
`pip install -r requirements.txt`, as above.

Note: scikit-build is only required for convenient management of the Python
build environment and packaging. See https://redmine.gromacs.org/issues/2896

# pybind11

Python bindings are expressed in C++ using the
[pybind11](https://pybind11.readthedocs.io/en/stable/)
template library.
The pybind11 repository is mirrored in GROMACS project sources and
installed with GROMACS for convenience and reproducibility.

# Build and install

## Cross compiling

On some systems, GROMACS will have been built and installed for a different
architecture than the system on which the Python package will be compiled.
We need to use CMake Tool-chains to support cross-compiling for the target architecture.

Note: scikit-build can use CMake Toolchains to properly handle `pip` builds.

## Offline installation

The `pip install` options `--no-index` and `--find-links` allow for an offline stash of package archives so that
satisfying dependencies for a new virtualenv does not require network access or lengthy build times.

# Dependencies

## OS X
Some dependencies (notably, a Python installation itself) may require some fiddling
with the XCode SDK.
https://developer.apple.com/documentation/xcode_release_notes/xcode_10_release_notes#3035624
