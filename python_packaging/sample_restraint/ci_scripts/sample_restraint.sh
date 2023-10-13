#!/usr/bin/env bash
set -ev

rm -rf build
mkdir build
pushd build
 cmake .. -DPYTHON_EXECUTABLE=$PYTHON
 make -j2 install
 make -j2 test
 $PYTHON -c "import myplugin"
popd
pushd tests
 $PYTHON -m pytest
 mpiexec -n 2 $PYTHON -m mpi4py -m pytest --log-cli-level=DEBUG -s --verbose
popd
