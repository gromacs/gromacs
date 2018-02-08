#!/bin/bash

#LOC="/opt/mackerell/apps/gromacs/drude"
LOC="$HOME/software/gromacs/git-drude_test5"
#CMAKE="/home/jalemkul/software/cmake-2.8.12/bin/cmake"
CMAKE=`which cmake`

# get rid of previous install
#echo "Nuking old GROMACS installation in $LOC..."
#rm -rf $LOC

# configure with fully static linking to all external libraries
CFLAGS="-static -static-libgcc" CXXFLAGS="-static -static-libstdc++" $CMAKE ../ -DCMAKE_PREFIX_PATH=/opt/fftw/3.3.4 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_INSTALL_PREFIX=$LOC -DGMX_GPU=OFF -DBUILD_SHARED_LIBS=OFF -DGMX_BUILD_SHARED_EXE=OFF 

#-DCMAKE_BUILD_TYPE=Debug

# compile
make -j 4

# install
make install

# clean
rm -rf `ls | grep -v install_gmx.sh`

exit;
