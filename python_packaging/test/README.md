# Tests for gmxapi Python packages distributed with GROMACS

Potentially longer-running than the unit tests in the Python package directory
and broader in scope. Intended to include more thorough and complete integration
testing than in the `acceptance` subdirectory.

## Requirements

Python tests use the `unittest` standard library module and the `unittest.mock`
submodule, included in Python 3.3+.

The additional `pytest` package allows tests to be written more easily and
concisely, with easier test fixtures (through decorators), log handling, and
other output handling.

## Files

Python files beginning with `test_` are collected by the Python testing
framework during automatic test discovery.

`conftest.py` and `pytest.ini` provide configuration for `pytest`.

Tests are organized according to the functional requirements
documented in `roadmap.rst` in the `python_packaging` directory.

## Usage

For basic tests, install the Python package(s) (presumably into a virtualenv),
then use the `pytest` executable to run these tests against the installed
package(s).

`pytest $LOCAL_REPO_DIR/python_packaging/test`

where `$LOCAL_REPO_DIR` is the path to the local copy of the GROMACS source repository.

For multi-process tests, run with an MPI execution wrapper and the `mpi4py` module.

`mpiexec -n 2 python -m mpi4py -m pytest $LOCAL_REPO_DIR/python_packaging/test`
