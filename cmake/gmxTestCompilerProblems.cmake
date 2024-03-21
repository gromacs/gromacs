#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

include(CheckCXXSourceCompiles)

# Macro that runs through a number of tests for buggy compiler
# versions, or other potential problems.
macro(gmx_test_compiler_problems)

    # Warn if C and C++ compilers do not match
    if(NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_CXX_COMPILER_ID)
        message(WARNING "The ids of the C and C++ compilers do not match (${CMAKE_C_COMPILER_ID} and ${CMAKE_CXX_COMPILER_ID}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()
    if(NOT CMAKE_C_COMPILER_VERSION STREQUAL CMAKE_CXX_COMPILER_VERSION)
        message(WARNING "The versions of the C and C++ compilers do not match (${CMAKE_C_COMPILER_VERSION} and ${CMAKE_CXX_COMPILER_VERSION}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()

    # Error if compiler doesn't support required C++17 features.
    # cmake feature detection is currently inconsistent: gitlab.kitware.com/cmake/cmake/issues/18869
    # We might want to switch to using feature test macros some time.
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
            set(cxx_required_version "GCC version 9")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19.15)
            set(cxx_required_version "Visual Studio 2017")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5)
            set(cxx_required_version "Clang 5")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        # All versions of IntelLLVM (a.k.a. DPCPP) compiler so far support C++17
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
        # All versions of NVC++ compiler so far support C++17
        # but only due to other bugs related to OMP 24.1+ version works
        # by disabling some OMP regions for it.
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 24.1)
            set(cxx_required_version "NVHPC SDK 24.1")
        endif()
    else()
        message(WARNING "You are using an unsupported compiler. Please make sure it fully supports C++17.")
    endif()
    if (cxx_required_version)
        message(FATAL_ERROR "${cxx_required_version} or later required. "
                            "Earlier versions don't have full C++17 support.")
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        message(WARNING "The Intel classic compiler is no longer supported. It may pass the tests, but is not tested by the GROMACS developers. Use the clang-based compiler from oneAPI, or gcc")
    endif()
    # Intel LLVM 2021.2 defaults to no-finite-math which isn't OK for GROMACS and its dependencies (muParser and GTest).
    # This is why we set the flags globally via CMAKE_CXX_FLAGS
    if(GMX_INTEL_LLVM AND GMX_INTEL_LLVM_VERSION GREATER_EQUAL 2021020)
        if (NOT WIN32)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-finite-math-only")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
        endif()
    endif()


    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
        check_cxx_source_compiles(
"// Test in-class array initalizers used with constructor initializer lists
struct TestStruct
{
    float a[3][3] = {{0}}; // in-class initializer
    float b; // not initialized until constructor initializer list
    TestStruct();
};
TestStruct::TestStruct() : b(0) {}
}" XLC_COMPILES_CORRECTLY)
        if (NOT XLC_COMPILES_CORRECTLY)
            message(FATAL_ERROR "No known version of xlC can compile the normal C++11 code in GROMACS, highest version checked is 16.1.0")
        endif()
    endif()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
        message(WARNING "Currently tested PGI compiler versions (up to 15.7) generate binaries that do not pass all regression test, and the generated binaries are significantly slower than with GCC or Clang. For now we do not recommend PGI beyond development testing - make sure to run the regressiontests.")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        check_cxx_source_compiles(
"#if !defined(__cray__) || !__cray__
#error Not Cray
#endif
int main(void){return 0;}" HAVE_CRAY_MACRO)
        # Later we will check if we can turn off GMX_THREAD_MPI but this check needs to run first
        if (HAVE_CRAY_MACRO AND GMX_OPENMP AND GMX_THREAD_MPI AND NOT GMX_MPI)
            message(WARNING "Clang-based Cray compiler works poorly with OpenMP and threadMPI. Consider enabling libMPI with -DGMX_MPI=ON.")
        endif()
    endif()


endmacro(gmx_test_compiler_problems)
