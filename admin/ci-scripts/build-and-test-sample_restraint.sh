#!/usr/bin/env bash
#
# Build and test the sample_restraint package distributed with GROMACS.
#
# This script is intended to support automated GROMACS testing infrastructure,
# and may be removed without notice.
#
# WARNING: This script assumes OpenMPI mpiexec. Syntax for launch wrappers from
# other implementations will need different syntax, and we should get a
# MPIRUNNER from the environment, or something.

# Make sure the script errors if any commands error.
set -e

# In the build-and-test-py-gmxapi.sh script, we use the currently canonical mechanism for hinting
# gromacs client build systems, but users have become accustomed to simply sourcing the GMXRC, so
# we test that alternative when building the gmxapi Python package here.
. $GROMACS_ROOT/bin/GMXRC

# Build and install gmxapi python package from local source.
# Note that tool chain may be provided differently across GROMACS versions.
python -m pip install --upgrade scikit-build-core[pyproject]
python -m pip install --no-build-isolation --no-cache-dir --no-deps --no-index python_packaging/gmxapi
# For more output, use the following instead
#python -m pip install --verbose --no-build-isolation --no-cache-dir --no-deps --no-index python_packaging/gmxapi --config-settings=cmake.verbose=true --config-settings=logging.level=DEBUG

pushd python_packaging/sample_restraint
  rm -rf build
  mkdir build
  pushd build
    # TODO: Update with respect to https://gitlab.com/gromacs/gromacs/-/issues/3133
    cmake .. \
      -C $GROMACS_ROOT/share/cmake/gromacs${GMX_SUFFIX}/gromacs-hints${GMX_SUFFIX}.cmake \
      -DPYTHON_EXECUTABLE=`which python` \
      -DGMXAPI_EXTENSION_DOWNLOAD_PYBIND=ON
    make -j4 tests
    make test
    #TODO: Can we get ctest JUnitXML output here?

    make install
  popd

  python -m pytest $PWD/tests --junitxml=$PLUGIN_TEST_XML --threads=${KUBERNETES_CPU_REQUEST}

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
        -n $((KUBERNETES_CPU_REQUEST/2)) ${PROGRAM[@]} --junitxml=$PLUGIN_MPI_TEST_XML : \
        -n $((KUBERNETES_CPU_REQUEST/2)) ${PROGRAM[@]}
  fi
popd
