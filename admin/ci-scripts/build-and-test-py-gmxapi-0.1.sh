#!/usr/bin/env bash
#
# Build, install, and test the gmxapi 0.1.x Python package distributed with
# GROMACS 2020.
#
# It is assumed that this script is run in an activated Python venv with the
# gmxapi dependencies already installed, with `python` resolvable by the shell
# to the appropriate Python interpreter.
#
# This script is intended to support automated GROMACS testing infrastructure,
# and may be removed without notice.
#
# WARNING: This script assumes OpenMPI mpiexec. Syntax for launch wrappers from
# other implementations will need different syntax, and we should get a
# MPIRUNNER from the environment, or something.

# Make sure the script errors if any commands error.
set -ev

# Create "sdist" source distribution archive.
pushd python_packaging/src
  GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
      python setup.py sdist
  # TODO: Identify SDIST

  # Build and install from sdist.
  GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
      python -m pip install \
          --no-cache-dir \
          --no-deps \
          --no-index \
          --no-build-isolation \
          dist/gmxapi*
  # TODO: Build and install from $SDIST instead of wildcard.
popd

# Run Python unit tests.
python -m pytest $PWD/python_packaging/src/test --junitxml=$PY_UNIT_TEST_XML
# TODO: enable MPI tests
#if [ -x `which mpiexec` ]; then
#    PYTHONDONTWRITEBYTECODE=1 \
#    mpiexec --allow-run-as-root \
#      --mca opal_warn_on_missing_libcuda 0 \
#      --mca orte_base_help_aggregate 0 \
#      -n 2 \
#      `which python` -m pytest \
#        -p no:cacheprovider \
#        $PWD/python_packaging/src/test \
#        --junitxml=$PY_MPI_UNIT_TEST_XML
#fi

# Run Python acceptance tests.
python -m pytest $PWD/python_packaging/test --junitxml=$PY_ACCEPTANCE_TEST_XML
# TODO: enable MPI tests
#if [ -x `which mpiexec` ]; then
#    PYTHONDONTWRITEBYTECODE=1 \
#    mpiexec --allow-run-as-root \
#      --mca opal_warn_on_missing_libcuda 0 \
#      --mca orte_base_help_aggregate 0 \
#      -n 2 \
#      `which python` -m pytest \
#        -p no:cacheprovider \
#        $PWD/python_packaging/test \
#        --junitxml=$PY_MPI_ACCEPTANCE_TEST_XML
#fi
