#!/bin/bash

set -ev

SCRIPT=$PWD/scripted_gmx_docker_builds.py

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

tag="gromacs/cmake-3.9.6-gcc-5-cuda-9.0-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.9.6 --gcc 5 --cuda 9.0 --ubuntu 16.04 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.9.6-gcc-6-cuda-10.1-nvidiaopencl-clfft-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.9.6 --gcc 6 --cuda 10.1 --opencl --clfft --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.9.6-gcc-7-amdopencl-clfft-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.9.6 --gcc 7 --opencl amd --clfft --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.15.7-gcc-8-cuda-10.1-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --gcc 8 --cuda 10.1 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.9.6-gcc-9-cuda-10.0-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.9.6 --gcc 9 --cuda 10.0 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.11.4-llvm-8-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.11.4 --llvm 8 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-8-tsan:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 8 --tsan | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-8-cuda-10.1-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 8 --cuda 10.1 --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.15.7-llvm-8-intelopencl-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --cmake 3.15.7 --llvm 8 --opencl intel --mpi openmpi | docker build -t $tag -

tag="gromacs/cmake-3.9.6-llvm-3.6-amdopencl-openmpi:2020"
tags[${#tags[@]}]=$tag
python3 $SCRIPT --ubuntu 16.04 --cmake 3.9.6 --llvm 3.6 --opencl amd --mpi openmpi | docker build -t $tag -

tag=gromacs/ci-docs-llvm:2020
tags[${#tags[@]}]=$tag
python3 $SCRIPT --llvm --doxygen | docker build -t $tag -

tag=gromacs/ci-docs-gcc:2020
tags[${#tags[@]}]=$tag
python3 $SCRIPT --gcc --doxygen | docker build -t $tag -

docker login
for tag in "${tags[@]}"; do
  echo "Pushing $tag"
  #docker push $tag
done
