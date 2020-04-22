#!/usr/bin/env bash
#
# Build and test the sample_restraint package distributed with GROMACS 2021.
#
# This script is intended to support automated GROMACS testing infrastructure,
# and may be removed without notice.
#
# WARNING: This script assumes OpenMPI mpiexec. Syntax for launch wrappers from
# other implementations will need different syntax, and we should get a
# MPIRUNNER from the environment, or something.

# Make sure the script errors if any commands error.
set -ev

source $VENVPATH/bin/activate

# Create "sdist" source distribution archive.
pushd python_packaging/src
  # TODO: Remove extraneous environment variable with resolution of #3273
  # Ref: https://redmine.gromacs.org/issues/3273
  GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
      python setup.py sdist
  # TODO: Identify SDIST

  # Build and install from sdist.
  # Note that tool chain may be provided differently in GROMACS 2020 and 2021.
  # TODO: How should we distinguish the target branch?
#  if [ "2020" -eq "$GROMACS_VERSION" ]; then
  if [ 1 ]; then
      GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
          python -m pip install \
              --no-cache-dir \
              --no-deps \
              --no-index \
              --no-build-isolation \
              dist/gmxapi*
  # TODO: Build and install from $SDIST instead of wildcard.
  # TODO: How should we distinguish the target branch?
#    elif [ "2021" -eq "$GROMACS_VERSION" ]; then
#      GMXTOOLCHAINDIR=$INSTALL_DIR/share/cmake/gromacs \
#          python -m pip install \
#              --no-cache-dir \
#              --no-deps \
#              --no-index \
#              --no-build-isolation \
#              dist/gmxapi*
  # TODO: Build and install from $SDIST instead of wildcard.
  else
      echo "Logic error in GROMACS version handling."
      exit 1
  fi
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
# TODO: enable MPI tests
#  if [ -x `which mpiexec` ]; then
#      PYTHONDONTWRITEBYTECODE=1 \
#      mpiexec --allow-run-as-root \
#        --mca opal_warn_on_missing_libcuda 0 \
#        --mca orte_base_help_aggregate 0 \
#        -n 2 \
#        `which python` -m pytest \
#          -p no:cacheprovider \
#          $PWD/tests \
#          --junitxml=$PLUGIN_MPI_TEST_XML
#  fi
popd
