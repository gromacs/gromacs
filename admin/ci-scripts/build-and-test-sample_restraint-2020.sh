#!/usr/bin/env bash
#
# Build and test the sample_restraint package distributed with GROMACS 2020.
#
# This script is intended to support automated GROMACS testing infrastructure,
# and may be removed without notice.
#
# WARNING: This script assumes OpenMPI mpiexec. Syntax for launch wrappers from
# other implementations will need different syntax, and we should get a
# MPIRUNNER from the environment, or something.

# Make sure the script errors if any commands error.
set -e

# TODO: We should be able to just use a $GMXAPI_0_1_SDIST or venv artifact...
pushd python_packaging/src
  GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
    python -m pip install \
        --no-cache-dir \
        --no-deps \
        --no-index \
        --no-build-isolation \
        .
popd

. $INSTALL_DIR/bin/GMXRC
pushd python_packaging/sample_restraint
  mkdir build
  pushd build
    # TODO: Update with respect to https://redmine.gromacs.org/issues/3133
    cmake .. \
             -DDOWNLOAD_GOOGLETEST=ON \
             -DGMXAPI_EXTENSION_DOWNLOAD_PYBIND=ON
    make

    make test
    #TODO: Can we get ctest JUnitXML output here?

    make install
  popd

  python -m pytest $PWD/tests --junitxml=$PLUGIN_TEST_XML

  # Note: Multiple pytest processes getting --junitxml output file argument
  # may cause problems, so we set the option on only one of the launched processes.
  # See also Multiple Instruction Multiple Data Model for OpenMPI mpirun:
  # https://www.open-mpi.org/doc/v3.0/man1/mpiexec.1.php
  PROGRAM=(`which python` -m mpi4py -m pytest \
          -p no:cacheprovider \
          $PWD/tests)
  # shellcheck disable=SC2068
  if [ -x `which mpiexec` ]; then
      PYTHONDONTWRITEBYTECODE=1 \
      mpiexec --allow-run-as-root \
        --mca opal_warn_on_missing_libcuda 0 \
        --mca orte_base_help_aggregate 0 \
        -n 1 ${PROGRAM[@]} --junitxml=$PLUGIN_MPI_TEST_XML : \
        -n 1 ${PROGRAM[@]}
  fi
popd
