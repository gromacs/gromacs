#!/bin/bash

set -ev

SCRIPT=$PWD/scripted_gmx_docker_builds.py

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

args[${#args[@]}]="--gcc 10 --clfft --mpi openmpi --rocm"
args[${#args[@]}]="--gcc 9 --clfft --mpi openmpi --rocm"
args[${#args[@]}]="--gcc 10 --cuda 11.2.2 --clfft --mpi openmpi"
args[${#args[@]}]="--gcc 7 --cuda 11.0 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 11 --tsan"
args[${#args[@]}]="--llvm 8 --cuda 11.0 --clfft --mpi openmpi"
args[${#args[@]}]="--llvm 9 --clfft --mpi openmpi --rocm"
args[${#args[@]}]="--oneapi 2021.1.1"
args[${#args[@]}]="--oneapi 2021.2.0 --intel-compute-runtime 21.21.19914"
args[${#args[@]}]="--llvm --doxygen --mpi openmpi --venvs 3.7.7"
args[${#args[@]}]="--llvm 11 --cuda 11.2.2 --hipsycl 0bf6420aab18 --rocm 4.2"

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
