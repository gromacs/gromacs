cmake_minimum_required(VERSION 3.9.6)

# This needs to be set before project() in order to pick up toolchain files
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# OS X deployment target should be >=10.9 for modern C++ compatibility.
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.9 CACHE STRING
    "OS X deployment target below 10.9 does not use C++11 standard library"
    FORCE)
set(CMAKE_OSX_ARCHITECTURES x86_64 CACHE STRING
    "Python package should only be built in 64-bit architecture"
    FORCE)

project(gmxapi VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0054 NEW)
# honor the language standard settings for try_compile()
cmake_policy(SET CMP0067 NEW)

# It may help to disable extensions to encourage consistent standard library implementations,
# but I don't know whether this would disable necessary attributes for symbol visibility
# or unused variable annotation, etc. TBD. Also, may be incompatible with pybind config (also TBD).
# See https://crascit.com/2015/03/28/enabling-cxx11-in-cmake/
#set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(external/pybind)

pybind11_add_module(_gmxapi gmxapi/_gmxapi.cpp)

install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)

if(SKBUILD)
    install(TARGETS _gmxapi LIBRARY DESTINATION gmxapi)
else()
    message(WARNING "CMake should be run by scikit-build through setup.py to get pure Py sources.")
endif()
