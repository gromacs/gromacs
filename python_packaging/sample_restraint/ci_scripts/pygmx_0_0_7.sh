#!/usr/bin/env bash
set -ev

pushd $HOME
 [ -d gmxapi ] || git clone --depth=1 --no-single-branch https://github.com/kassonlab/gmxapi.git
 pushd gmxapi
  git checkout release-0_0_7
  rm -rf build
  mkdir -p build
  pushd build
   cmake .. -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_C_COMPILER=$CC -DPYTHON_EXECUTABLE=$PYTHON
   make -j2 install
  popd
 popd
 mpiexec -n 2 $PYTHON -m mpi4py -m pytest --log-cli-level=WARN --pyargs gmx -s
# mpiexec -n 2 $PYTHON -m mpi4py -m pytest --log-cli-level=DEBUG --pyargs gmx -s --verbose
 ccache -s
popd
