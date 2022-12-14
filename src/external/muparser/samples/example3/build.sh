#!/bin/bash -x

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MP_SOURCES=${CWD}/../../
MP_BUILD=${CWD}/muparser-build
MP_INSTALL=${CWD}/muparser-install
EX3_BUILD_TREE=${CWD}/example3-using-buildtree
EX3_INSTALL_TREE=${CWD}/example3-using-installtree

# Build muparser and install it
cmake -H${MP_SOURCES} -B${MP_BUILD} -DCMAKE_INSTALL_PREFIX=${MP_INSTALL}
cmake --build ${MP_BUILD} --target install

# Build the example using muparser build tree
cmake -H${CWD} -B${EX3_BUILD_TREE} -DCMAKE_PREFIX_PATH=${MP_BUILD} 
cmake --build ${EX3_BUILD_TREE} --target all
cmake --build ${EX3_BUILD_TREE} --target test

# Build the example using muparser install tree
cmake -H${CWD} -B${EX3_INSTALL_TREE} -DCMAKE_PREFIX_PATH=${MP_INSTALL} 
cmake --build ${EX3_INSTALL_TREE} --target all
cmake --build ${EX3_INSTALL_TREE} --target test

