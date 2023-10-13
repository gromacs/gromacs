#!/usr/bin/env bash
set -ev

export GMX_DOUBLE=OFF
export GMX_MPI=OFF
export GMX_THREAD_MPI=ON

export GMX_SRC_DIR=gromacs-2019

export CCACHE_DIR=$HOME/.ccache_gmxapi
ccache -s

pushd $HOME
 [ -d gromacs-gmxapi ] || \
    git clone \
        --depth=1 \
        --no-single-branch \
        https://github.com/gromacs/gromacs.git \
        ${GMX_SRC_DIR}
 pushd ${GMX_SRC_DIR}
  git branch -a
  git checkout release-2019
  pwd
  rm -rf build
  mkdir build
  pushd build
   cmake -DCMAKE_CXX_COMPILER=$CXX \
         -DGMX_ENABLE_CCACHE=ON \
         -DCMAKE_C_COMPILER=$CC \
         -DGMX_DOUBLE=$GMX_DOUBLE \
         -DGMX_MPI=$GMX_MPI \
         -DGMX_THREAD_MPI=$GMX_THREAD_MPI \
         -DGMXAPI=ON \
         -DCMAKE_INSTALL_PREFIX=$HOME/install/gromacs_2019 \
         ..
   make -j2 install
  popd
 popd
popd
ccache -s
