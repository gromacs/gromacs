#!/usr/bin/env bash
#
# Build and test the sample_restraint package distributed with GROMACS 2022.
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
  # Build and install gmxapi python package from local source.
  # Note that tool chain may be provided differently across GROMACS versions.
  if [ "2022" -eq "$GROMACS_MAJOR_VERSION" ]; then
      GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
          python -m pip install \
              --no-cache-dir \
              .
  # TODO: Get sdist artifact for other gmxapi versions.
#    elif [ "2021" -eq "$GROMACS_MAJOR_VERSION" ]; then
#      GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
#          python -m pip install \
#              --no-cache-dir \
#              --no-deps \
#              --no-index \
#              dist/gmxapi*
  else
      echo "Logic error in GROMACS version handling."
      exit 1
  fi
popd

. $INSTALL_DIR/bin/GMXRC
pushd python_packaging/sample_restraint
  mkdir build
  pushd build
    # TODO: Update with respect to https://gitlab.com/gromacs/gromacs/-/issues/3133
    cmake .. \
      -DPYTHON_EXECUTABLE=`which python` \
      -DDOWNLOAD_GOOGLETEST=ON \
      -DGMXAPI_EXTENSION_DOWNLOAD_PYBIND=ON
    make -j4 tests
    make test
    #TODO: Can we get ctest JUnitXML output here?

    make install
  popd

  python -m pytest $PWD/tests --junitxml=$PLUGIN_TEST_XML --threads=2

  # Note: Multiple pytest processes getting --junitxml output file argument
  # may cause problems, so we set the option on only one of the launched processes.
  # See also Multiple Instruction Multiple Data Model for OpenMPI mpirun:
  # https://www.open-mpi.org/doc/v3.0/man1/mpiexec.1.php
  PROGRAM=(`which python` -m mpi4py -m pytest \
          -p no:cacheprovider \
          $PWD/tests \
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
popd
