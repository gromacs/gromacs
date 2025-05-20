#!/usr/bin/env bash

set -e

SCRIPT=$PWD/scripted_gmx_docker_builds.py
PYTHON=${PYTHON:-$(which python3)}

# Note: All official GROMACS CI images are built
# with openmpi on. That reduces the total number of
# images needed, because the same one can test library,
# thread and no MPI configurations.

# These rows correspond to the containers that are currently used in
# GROMACS CI testing for this source-code branch. Release branches
# seldom ever change these, because only bug fixes are accepted and
# minimum support levels seldom change.
#
# In main branch, these are changed often as new minimum support
# levels are adopted and new capabilities of dependencies are
# utilized.
#
# This script outputs bash commands that are useful to make
# all GROMACS CI containers. In the usual case of not re-making
# all containers, it can be useful to use grep to select the
# relevant commands and pipe that to bash

args[${#args[@]}]="--ubuntu 22.04 --gcc 12 --clfft --mpi openmpi --rocm 5.4.1 --hdf5"
args[${#args[@]}]="--ubuntu 22.04 --gcc 13 --cuda 12.5.1 --clfft --mpi openmpi --nvhpcsdk 24.7"
args[${#args[@]}]="--ubuntu 22.04 --gcc 12 --cuda 12.1.0 --clfft --mpi openmpi --heffte v2.4.0 --libtorch"
args[${#args[@]}]="--ubuntu 24.04 --gcc 14 --mpi openmpi --cp2k 2024.2"
args[${#args[@]}]="--ubuntu 24.04 --gcc 11 --mpi openmpi --cp2k 9.1"
args[${#args[@]}]="--ubuntu 22.04 --llvm 18 --cuda 12.1.0"
args[${#args[@]}]="--ubuntu 24.04 --llvm 18 --tsan"
args[${#args[@]}]="--ubuntu 22.04 --llvm 14 --cuda 12.1.0 --clfft --mpi openmpi"
args[${#args[@]}]="--ubuntu 24.04 --llvm 19 --mpi openmpi --hdf5"
args[${#args[@]}]="--ubuntu 22.04 --oneapi 2025.1 --intel-compute-runtime"
args[${#args[@]}]="--ubuntu 22.04 --oneapi 2025.0 --rocm 6.1.3 --cuda 12.0.1 --oneapi-plugin-amd --oneapi-plugin-nvidia"
args[${#args[@]}]="--ubuntu 24.04 --llvm 19 --doxygen --mpi openmpi --venvs 3.9.13 3.12.5"
args[${#args[@]}]="--ubuntu 24.04 --llvm 18 --cuda 12.6.3 --adaptivecpp 24.10.0 --rocm 6.3.1 --mpi mpich"
args[${#args[@]}]="--ubuntu 22.04 --adaptivecpp 23.10.0 --rocm 5.7.1"
args[${#args[@]}]="--ubuntu 24.04 --rocm 6.2.2 --mpi openmpi --plumed --heffte v2.4.0"

echo
echo "Consider pulling the following images for build layer cache."
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  echo "docker pull $($PYTHON -m utility $arg_string)"
done
echo
echo

echo "To build with cache hints:"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$($PYTHON -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  # shellcheck disable=SC2086
  echo "$PYTHON $SCRIPT $arg_string | docker build -t $tag --cache-from $tag -"
done
unset tags
echo
echo

echo "To build without extra cache hints:"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$($PYTHON -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  # shellcheck disable=SC2086
  echo "$PYTHON $SCRIPT $arg_string | docker build -t $tag -"
done
unset tags
echo
echo

echo "Run the following to upload the updated images."
echo "docker login registry.gitlab.com -u <token name> -p <hash>"
echo
for arg_string in "${args[@]}"; do
  # shellcheck disable=SC2086
  tag=$($PYTHON -m utility $arg_string)
  tags[${#tags[@]}]=$tag
  echo "docker push $tag"
done

# Check whether built images are used and whether used images are built.
# Note: Let used images with a ':latest' tag match images without explicit tags.
for tag in "${tags[@]}"; do
  image=$(basename $tag)
  # Checking whether $image is used.
  grep -qR -e "${image}\(:latest\)*$" ../gitlab-ci || echo Warning: Image $image appears unused.
done
list=$(grep -R 'image: ' ../gitlab-ci/ |awk '{print $3}' |sort -u)
$PYTHON << EOF
from os.path import basename
built="""${tags[@]}"""
built=set(basename(image.rstrip()) for image in built.split())
in_use="""${list[@]}"""
in_use=[basename(image.rstrip()) for image in in_use.split()]
for tag in in_use:
  if tag.endswith(':latest'):
    tag = tag.split(':')[0]
  if not tag in built:
    print(f'Warning: {tag} is not being built.')
EOF
