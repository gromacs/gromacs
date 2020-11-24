#!/bin/bash

set -ev

SCRIPT=$PWD/scripted_gmx_docker_builds.py

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

args[${#args[@]}]="--gcc 8 --cuda 11.0 --clfft --mpi openmpi"
args[${#args[@]}]="--gcc 7 --clfft --mpi openmpi --ubuntu 18.04"
args[${#args[@]}]="--llvm 8 --tsan"
args[${#args[@]}]="--llvm 8 --cuda 10.0 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 8 --cuda 10.1 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 8 --cuda 11.0 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 9 --clfft --mpi openmpi --ubuntu 18.04"
args[${#args[@]}]="--oneapi 2021.1-beta09"
args[${#args[@]}]="--llvm --doxygen"

echo "Building the following images."
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  python3 -m utility $arg_string
done
echo

for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$(python3 -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  # shellcheck disable=SC2086
  python3 $SCRIPT $arg_string | docker build -t $tag -
done

echo "Run the following to upload the updated images."
echo "docker login registry.gitlab.com -u <token name> -p <hash>"
for tag in "${tags[@]}"; do
  echo "docker push $tag"
done
