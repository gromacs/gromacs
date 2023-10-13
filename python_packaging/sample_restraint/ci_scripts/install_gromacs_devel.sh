#!/usr/bin/env bash
set -ev

export GMX_DOUBLE=OFF
export GMX_MPI=OFF
export GMX_THREAD_MPI=ON

export GMX_SRC_DIR=gromacs-kassonlab

export CCACHE_DIR=$HOME/.ccache_gmxapi
ccache -s

pushd $HOME
 [ -d "${GMX_SRC_DIR}" ] || \
     git clone \
         --depth=1 \
         --no-single-branch \
         https://github.com/kassonlab/gromacs-gmxapi.git \
         ${GMX_SRC_DIR}
 pushd ${GMX_SRC_DIR}
  git branch -a
  git checkout devel
  pwd
  rm -rf build
  mkdir build
  pushd build
   cmake -DGMX_BUILD_HELP=OFF \
         -DGMX_ENABLE_CCACHE=ON \
         -DCMAKE_CXX_COMPILER=$CXX \
         -DCMAKE_C_COMPILER=$CC \
         -DGMX_DOUBLE=$GMX_DOUBLE \
         -DGMX_MPI=$GMX_MPI \
         -DGMX_THREAD_MPI=$GMX_THREAD_MPI \
         -DCMAKE_INSTALL_PREFIX=$HOME/install/gromacs_devel \
         ..
   make -j2 install
  popd
 popd
popd
ccache -s
