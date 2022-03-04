#!/bin/bash

set -e

SCRIPT=$PWD/scripted_gmx_docker_builds.py

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

args[${#args[@]}]="--gcc 11 --clfft --mpi openmpi --rocm"
args[${#args[@]}]="--gcc 11 --cuda 11.4.1 --clfft --mpi openmpi --heffte v2.2.0"
args[${#args[@]}]="--gcc 7 --cuda 11.0 --clfft --mpi openmpi --heffte v2.2.0"
args[${#args[@]}]="--gcc 9 --mpi openmpi --cp2k 8.2"
args[${#args[@]}]="--gcc 9 --mpi openmpi --cp2k 9.1"
args[${#args[@]}]="--llvm 11 --cuda 11.4.1"
args[${#args[@]}]="--llvm 11 --tsan"
args[${#args[@]}]="--llvm 8 --cuda 11.0 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 13 --clfft --mpi openmpi --rocm"
args[${#args[@]}]="--oneapi 2021.4.0"
args[${#args[@]}]="--oneapi 2021.4.0 --intel-compute-runtime"
args[${#args[@]}]="--llvm --doxygen --mpi openmpi --venvs 3.7.7"
args[${#args[@]}]="--llvm 12 --cuda 11.4.3 --hipsycl c1246fd --rocm 5.0"

echo
echo "Consider pulling the following images for build layer cache."
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  echo "docker pull $(python3 -m utility $arg_string)"
done
echo
echo

echo "To build with cache hints:"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$(python3 -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  # shellcheck disable=SC2086
  echo "$(which python3) $SCRIPT $arg_string | docker build -t $tag --cache-from $tag -"
done
unset tags
echo
echo

echo "To build without extra cache hints:"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$(python3 -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  # shellcheck disable=SC2086
  echo "$(which python3) $SCRIPT $arg_string | docker build -t $tag -"
done
unset tags
echo
echo

echo "Run the following to upload the updated images."
echo "docker login registry.gitlab.com -u <token name> -p <hash>"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$(python3 -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  echo "docker push $tag"
done
