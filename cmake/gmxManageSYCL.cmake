#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2020,2021, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# OpenCL required version: 1.2 or newer
set(REQUIRED_SYCL_MIN_VERSION_MAJOR 1)
set(REQUIRED_SYCL_MIN_VERSION_MINOR 2)
set(REQUIRED_SYCL_MIN_VERSION ${REQUIRED_SYCL_MIN_VERSION_MAJOR}.${REQUIRED_SYCL_MIN_VERSION_MINOR})

set(GMX_GPU_SYCL ON)

# CMake issue tracking the efforts to make a universal upstream module:
# https://gitlab.kitware.com/cmake/cmake/-/issues/21711

option(GMX_SYCL_HIPSYCL "Use hipSYCL instead of Intel/Clang for SYCL compilation" OFF)

if(GMX_DOUBLE)
    message(FATAL_ERROR "SYCL acceleration is not available in double precision")
endif()

include(gmxFindFlagsForSource)

if(GMX_SYCL_HIPSYCL)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(FATAL_ERROR "hipSYCL can only be built with Clang++ compiler")
    endif()
    set(HIPSYCL_CLANG "${CMAKE_CXX_COMPILER}")
    # -Wno-unknown-cuda-version because Clang-11 complains about CUDA 11.0-11.2, despite working fine with them.
    # -Wno-unknown-attributes because hipSYCL does not support reqd_sub_group_size (because it can only do some sub group sizes).
    # --hipsycl-explicit-multipass is needed when building for both CUDA and HIP.
    set(HIPSYCL_SYCLCC_EXTRA_ARGS "-Wno-unknown-cuda-version -Wno-unknown-attributes --hipsycl-explicit-multipass")
    find_package(hipsycl REQUIRED)
else()
    if(CMAKE_CXX_COMPILER MATCHES "dpcpp")
        # At least Intel dpcpp defaults to having SYCL enabled for all code. This leads to two problems:
        #
        # 1. Compiles take ~3x longer, since every file has to be compiled for multiple targets.
        # 2. We get a ton of warnings for the device-specific pass when the compiler sees our SIMD code.
        #
        # To avoid this, we attempt to find a flag to disable SYCL for non-SYCL files. Unfortunately,
        # when using gmx_find_flag_for_source() that includes calling check_cxx_compiler_flag(),
        # this in turn exposes a bug in dpcpp, where an object file compiles with -fno-sycl leads to
        # a failed link stage (when the same flag is not used). Since none of this is critical, we handle
        # it by merely checking if it works to compile a source fils with this flag, and choking if SYCL
        # is still enabled.
    
        if(NOT CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY)
            message(STATUS "Checking for flags to disable SYCL")
        endif()
    
        gmx_check_source_compiles_with_flags(
            "int main() { return 0; }"
            "-fno-sycl"
            "CXX"
            DISABLE_SYCL_CXX_FLAGS_RESULT)
    
        if(DISABLE_SYCL_CXX_FLAGS_RESULT)
            set(DISABLE_SYCL_CXX_FLAGS "-fno-sycl")
        endif()
        if(NOT CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY)
            if(DISABLE_SYCL_CXX_FLAGS_RESULT)
                message(STATUS "Checking for flags to disable SYCL - -fno-sycl")
            else()
                message(WARNING "Cannot find flags to disable SYCL for non-SYCL hardware-specific C++ code. Expect many warnings, but they are likely benign.")
            endif()
            set(CHECK_DISABLE_SYCL_CXX_FLAGS_QUIETLY 1 CACHE INTERNAL "Keep quiet on future calls to detect no-SYCL flags" FORCE)
        endif()
    endif()
    
    # Find the flags to enable (or re-enable) SYCL with Intel extensions. In case we turned it off above,
    # it's important that we check the combination of both flags, to make sure the second one re-enables SYCL.
    if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
        message(STATUS "Checking for flags to enable SYCL")
    endif()
    gmx_find_flag_for_source(SYCL_CXX_FLAGS_RESULT
        "#include <CL/sycl.hpp>
         namespace sycl = cl::sycl;
         int main(){
             sycl::queue q(sycl::default_selector{});
             return 0;
         }
         " "CXX" DISABLE_SYCL_CXX_FLAGS SYCL_CXX_FLAGS "-fsycl -fsycl-device-code-split=per_kernel")
    
    if(NOT CHECK_SYCL_CXX_FLAGS_QUIETLY)
        if(SYCL_CXX_FLAGS_RESULT)
            message(STATUS "Checking for flags to enable SYCL - ${SYCL_CXX_FLAGS}")
        endif()
        set(CHECK_SYCL_CXX_FLAGS_QUIETLY 1 CACHE INTERNAL "Keep quiet on future calls to detect SYCL flags" FORCE)
    endif()
    
    if(NOT SYCL_CXX_FLAGS_RESULT)
        message(FATAL_ERROR "Cannot compile with SYCL Intel compiler. Try a different compiler or disable SYCL.")
    endif()

    # Add function wrapper similar to the one used by ComputeCPP and hipSYCL
    function(add_sycl_to_target)
        cmake_parse_arguments(
            PARSE_ARGV 0 # No positional arguments
            ARGS # Prefix for the resulting variables
            "" # No options
            "TARGET" # One-value keyword
            "SOURCES" # Multi-value keyword
        )
        set_source_files_properties(${ARGS_SOURCES} PROPERTIES COMPILE_FLAGS "${SYCL_CXX_FLAGS}")
        target_link_libraries(${ARGS_TARGET} PRIVATE ${SYCL_CXX_FLAGS})
    endfunction(add_sycl_to_target)
endif()
