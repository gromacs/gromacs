#!/bin/bash

set -ev

SCRIPT=$PWD/scripted_gmx_docker_builds.py

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

tag="gromacs/cmake-3.15.7-gcc-8-cuda-10.1-nvidiaopencl-clfft-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --gcc 8 --cuda 10.1 --opencl --clfft --mpi openmpi \
| docker build -t $tag -

tag="gromacs/cmake-3.13.0-gcc-7-amdopencl-clfft-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.13.0 --gcc 7 --opencl amd --clfft --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.13.0-llvm-8-tsan:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.13.0 --llvm 8 --tsan | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-8-cuda-10.0-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 8 --cuda 10.0 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-8-cuda-10.1-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 8 --cuda 10.1 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-9-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 9 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.13.0-llvm-9-intelopencl-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.13.0 --llvm 9 --opencl intel --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.13.0-llvm-9-amdopencl-openmpi:master"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.13.0 --llvm 9 --opencl amd --mpi openmpi | docker build -t $tag -

tag=gromacs/ci-docs-llvm:master
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.17.2 --llvm --doxygen | docker build -t $tag -

tag=gromacs/ci-docs-gcc:master
tags[${#tags[@]}]=$tag
python3 $SCRIPT --gcc --doxygen | docker build -t $tag -

echo "Run the following to upload the updated images."
echo docker login
for tag in "${tags[@]}"; do
  echo docker push $tag
done
