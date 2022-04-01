#!/usr/bin/env bash

# Configure and build the sample code for the Trajectory Analysis Framework.
# Input variables:
#   CMAKE: Executable (path) for the CMake binary. (E.g. `$(which cmake)`)
#   GROMACS_DIR: absolute path to a GROMACS installation
#   GMX_SUFFIX: (optional) suffix as used by template/CMakeLists.txt (E.g. "_mpi")
#   GMX_DOUBLE: (optional) requested floating point precision for template/CMakeLists.txt
#   TMP_BUILD: (optional) the name of a directory to create and use for the build

set -e

CLIENT_CMAKE_ARGS="-C ${GROMACS_DIR}/share/cmake/gromacs${GMX_SUFFIX}/gromacs-hints${GMX_SUFFIX}.cmake \
 -Werror=dev \
 -DGROMACS_DIR=${GROMACS_DIR}"

if [ "${GMX_SUFFIX}" ] ; then
  CLIENT_CMAKE_ARGS="${CLIENT_CMAKE_ARGS} -DGMX_SUFFIX=${GMX_SUFFIX}"
fi

if [ "${GMX_DOUBLE}" ]; then
    CLIENT_CMAKE_ARGS="${CLIENT_CMAKE_ARGS} -DGMX_DOUBLE=${GMX_DOUBLE}"
fi

BUILD=${TMP_BUILD:-"taf-example-build"}
mkdir -p ${BUILD}
pushd ${BUILD}
  ${CMAKE} \
    ${GROMACS_DIR}/share/gromacs/template \
    ${CLIENT_CMAKE_ARGS}
  ${CMAKE} --build . --target all -- VERBOSE=1
  ./template --version
popd
