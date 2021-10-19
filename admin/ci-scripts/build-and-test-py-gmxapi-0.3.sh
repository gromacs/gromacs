#!/usr/bin/env bash
#
# Build, install, and test the gmxapi 0.3 Python package developed with
# GROMACS 2022.
#
# This script assumes an activated Python venv with the
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
set -e

pushd python_packaging/src
  # Make sure to delete any accidentally lingering build artifacts.
  rm -rf build dist
  # Build and install the gmxapi Python package.
  # TODO(#3273): Reduce requirements for `setup.py` `sdist` command and provide build artifact.
  GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
      python -m pip install \
          --no-build-isolation \
          --no-cache-dir \
          --no-deps \
          --no-index \
          .
popd

# Run Python unit tests.
python -m pytest python_packaging/src/test --junitxml=$PY_UNIT_TEST_XML --threads=2

# Note: Multiple pytest processes getting --junitxml output file argument
# may cause problems, so we set the option on only one of the launched processes.
# See also Multiple Instruction Multiple Data Model for OpenMPI mpirun:
# https://www.open-mpi.org/doc/v3.0/man1/mpiexec.1.php
PROGRAM=(`which python` -m mpi4py -m pytest \
        -p no:cacheprovider \
        $PWD/python_packaging/src/test \
        --threads=1)
# shellcheck disable=SC2068
if [ -x `which mpiexec` ]; then
    PYTHONDONTWRITEBYTECODE=1 \
    mpiexec --allow-run-as-root \
      -x OMP_NUM_THREADS=1 \
      --mca opal_warn_on_missing_libcuda 0 \
      --mca orte_base_help_aggregate 0 \
      -n 1 ${PROGRAM[@]} --junitxml=$PLUGIN_MPI_TEST_XML : \
      -n 1 ${PROGRAM[@]}
fi
