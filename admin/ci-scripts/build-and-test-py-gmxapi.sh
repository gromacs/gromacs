#!/usr/bin/env bash
#
# Build, install, and test the gmxapi Python package.
#
# Inputs (environment variables):
#
#     GROMACS_ROOT: absolute path to GROMACS installation prefix.
#     GMX_SUFFIX: (optional) suffix as used by template/CMakeLists.txt (E.g. "_mpi")
#     PY_UNIT_TEST_XML: XML output file path for basic pytest suite
#     PY_MPI_UNIT_TEST_XML: XML output file path for MPI pytest suite
#
# This script assumes an activated Python venv with the
# gmxapi dependencies already installed, with `python` resolvable by the shell
# to the appropriate Python interpreter.
#
# This script is intended to support automated GROMACS testing infrastructure,
# and may be removed without notice.
#
# WARNING: This script assumes OpenMPI mpiexec. Launch wrappers from
# other implementations will need different syntax, and we should get a
# MPIRUNNER from the environment, or something.

# Make sure the script errors if any commands error.
set -e

pushd python_packaging/gmxapi
  # Make sure to delete any accidentally lingering build artifacts.
  rm -rf build dist
  # Build and install the gmxapi Python package.
  # Use the documented mechanism for getting GROMACS installation and build system hints.
  # See docs/gmxapi/userguide/install.rst and docs/release-notes/2022/major/portability.rst
  CMAKE_ARGS="-Dgmxapi_ROOT=$GROMACS_ROOT -C $GROMACS_ROOT/share/cmake/gromacs${GMX_SUFFIX}/gromacs-hints${GMX_SUFFIX}.cmake" \
      python -m pip install \
          --no-build-isolation \
          --no-cache-dir \
          --no-deps \
          --no-index \
          .
popd

# Run Python unit tests.
python -m pytest python_packaging/gmxapi/test --junitxml="$PY_UNIT_TEST_XML" --threads=${KUBERNETES_CPU_REQUEST}

# Note: Multiple pytest processes getting --junitxml output file argument
# may cause problems, so we set the option on only one of the launched processes.
# See also Multiple Instruction Multiple Data Model for OpenMPI mpirun:
# https://www.open-mpi.org/doc/v3.0/man1/mpiexec.1.php
PROGRAM=("$(which python)" -m mpi4py -m pytest \
        -p no:cacheprovider \
        "$PWD"/python_packaging/gmxapi/test \
        --threads=1)
if [ -x "$(which mpiexec)" ]; then
    PYTHONDONTWRITEBYTECODE=1 \
    mpiexec --allow-run-as-root \
      -x OMP_NUM_THREADS=1 \
      --mca opal_warn_on_missing_libcuda 0 \
      --mca orte_base_help_aggregate 0 \
      -n $((KUBERNETES_CPU_REQUEST/2)) "${PROGRAM[@]}" --junitxml="$PY_MPI_UNIT_TEST_XML" : \
      -n $((KUBERNETES_CPU_REQUEST/2)) "${PROGRAM[@]}"
fi
