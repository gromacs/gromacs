# Python package sources

This directory hosts subtrees for (separable) Python packages that are
maintained as part of the GROMACS repository.

## Repository organization

This directory contains subtrees for GROMACS Python packages, some supporting
Docker infrastructure, and CMake infrastructure for building the packages as
part of a GROMACS build (for testing purposes or when building complete
documentation).

**TODO: Consider long term homes of these directory contents.**

The `gmxapi` Python package might be more maintainable as a separate
repository (embedded in GROMACS via `git submodule` or `git subtree` for
convenience) to more closely match the `sdist` archives (see
https://pypi.org/project/gmxapi/)

The `sample_restraint` repository would be more helpful to researchers as a
[cookiecutter](https://cookiecutter.readthedocs.io/en/stable/). The
cookiecutter template project would necessarily be a separate repository, and
it would make sense to reference the cookiecutter repository (such as with
`git submodule`) rather than duplicating the subtree in this subdirectory of
the GROMACS repository.

## gmxapi

Python framework for gmxapi high-level interface.

Refer to ../docs/gmxapi/userguide/install.rst
(https://manual.gromacs.org/current/gmxapi/userguide/install.html)
for regular installation instructions.

Use `pytest` to run unit tests and integration tests.

    cd python_packaging/gmxapi
    pip install -r requirements.txt
    pytest test

To build the complete GROMACS documentation, or to run Python package tests
as part of the GROMACS `make check`, configure CMake with
`-DGMX_PYTHON_PACKAGE=ON` to do a partial installation into the GROMACS build
tree.

For additional discussion on packaging and distribution, see
https://gitlab.com/gromacs/gromacs/-/issues/2896

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

# Tests for gmxapi Python packages distributed with GROMACS

## Requirements

`gmxapi/requirements.txt` includes sufficient packages for complete testing and
documentation builds:

    pip install -r python_packaging/gmxapi/requirements.txt

## Files

Python files beginning with `test_` are collected by the Python testing
framework during automatic test discovery.

`conftest.py` and `pytest.ini` provide configuration for `pytest`.

## Usage

For basic tests, install the Python package(s) (presumably into a virtualenv),
then use the `pytest` executable to run these tests against the installed
package(s).

`pytest $LOCAL_REPO_DIR/python_packaging/gmxapi/test`

where `$LOCAL_REPO_DIR` is the path to the local copy of the GROMACS source repository.

For multi-process tests, run with an MPI execution wrapper and the `mpi4py` module.

`mpiexec -n 2 python -m mpi4py -m pytest -x $LOCAL_REPO_DIR/python_packaging/test`

The `-x` flag causes pytest to exit at the first test failure. This is important
to prevent hangs in an MPI context (in which a failure on one rank may prevent
it to reach a collective call entered on another rank).

## Controlling output

Refer to pytest documentation for command line options to control the type and detail of output.
Some high level overview and basic tasks are online at https://docs.pytest.org/en/3.9.3/usage.html
but you should run `pytest -h` in a terminal to get the complete set of available options
(in particular, note *log_cli* and *log_level*).
