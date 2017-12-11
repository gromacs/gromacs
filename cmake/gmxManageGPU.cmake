#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

# If the user did not set GMX_GPU we'll consider this option to be
# in "auto" mode meaning that we will:
# - search for CUDA and set GMX_GPU=ON we find it
# - check whether GPUs are present
# - if CUDA is not found but GPUs were detected issue a warning
if (NOT DEFINED GMX_GPU)
    set(GMX_GPU_AUTO TRUE CACHE INTERNAL "GPU acceleration will be selected automatically")
endif()
option(GMX_GPU "Enable GPU acceleration" OFF)

option(GMX_CLANG_CUDA "Use clang for CUDA" OFF)
if (GMX_CLANG_CUDA)
    # CUDA 7.0 or later required, override req. version
    set(REQUIRED_CUDA_VERSION 7.0)
endif()

if(GMX_GPU AND GMX_DOUBLE)
    message(FATAL_ERROR "GPU acceleration is not available in double precision!")
endif()
if(GMX_GPU_AUTO AND GMX_DOUBLE)
    message(WARNING "GPU acceleration is not available in double precision, disabled!")
    set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
    set_property(CACHE GMX_GPU_AUTO PROPERTY VALUE OFF)
endif()

# detect GPUs in the build host machine
if ((GMX_GPU OR GMX_GPU_AUTO) AND NOT GMX_GPU_DETECTION_DONE)
    include(gmxDetectGpu)
    gmx_detect_gpu()
endif()

# We need to call find_package even when we've already done the detection/setup
if(GMX_GPU OR GMX_GPU_AUTO)
    if(NOT GMX_GPU AND NOT GMX_DETECT_GPU_AVAILABLE)
        # Stay quiet when detection has occured and found no GPU.
        # Noise is acceptable when there is a GPU or the user required one.
        set(FIND_CUDA_QUIETLY QUIET)
    endif()

    # Cmake tries to use the static cuda runtime by default,
    # but this leads to unusable GPU builds on OS X.
    if(APPLE)
        set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE STRING "Use the static version of the CUDA runtime library if available")
    endif()

    find_package(CUDA ${REQUIRED_CUDA_VERSION} ${FIND_CUDA_QUIETLY})
endif()

# Depending on the current vale of GMX_GPU and GMX_GPU_AUTO:
# - OFF, FALSE: Will skip this detection/setup.
# - OFF, TRUE : Will keep GMX_GPU=OFF if no CUDA is detected, but will assemble
#               a warning message which will be issued at the end of the
#               configuration if GPU(s) were found in the build system.
# - ON , FALSE: The user requested GPU build and this requires CUDA, so we will
#               fail if it is not available.
# - ON , TRUE : Can't happen (GMX_GPU=ON can only be user-set at this point)
if((GMX_GPU OR GMX_GPU_AUTO) AND NOT GMX_GPU_DETECTION_DONE)
    if (EXISTS ${CUDA_TOOLKIT_ROOT_DIR})
        set(CUDA_FOUND TRUE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    else()
        set(CUDA_FOUND FALSE CACHE INTERNAL "Whether the CUDA toolkit was found" FORCE)
    endif()

    # assemble warning/error message
    if (GMX_DETECT_GPU_AVAILABLE)
        set(_msg "${GMX_DETECT_GPU_COUNT} NVIDIA GPU(s) found in the system")

        # append GPU names
        if (NOT GMX_DETECT_GPU_INFO STREQUAL "")
            set(_msg "${_msg}:")
            foreach(gpu ${GMX_DETECT_GPU_INFO})
                set(_msg "${_msg}
${gpu}")
            endforeach()
        endif()

        # TODO remove the second part of the message when we'll have compute
        # capability information from the detection.
        set(_msg "${_msg}
Compute capability information not available, consult the NVIDIA website:
https://developer.nvidia.com/cuda-gpus")
    endif()

        set(CUDA_NOTFOUND_MESSAGE "mdrun supports native GPU acceleration on NVIDIA hardware with compute capability >= ${REQUIRED_CUDA_COMPUTE_CAPABILITY} (Fermi or later). This requires the NVIDIA CUDA toolkit, which was not found. Its location can be hinted by setting the CUDA_TOOLKIT_ROOT_DIR CMake option (does not work as an environment variable). The typical location would be /usr/local/cuda[-version]. Note that CPU or GPU acceleration can be selected at runtime.

${_msg}")
        unset(_msg)

    if (NOT CUDA_FOUND)
        if (GMX_GPU_AUTO)
            # Disable GPU acceleration in auto mode
            message(STATUS "No compatible CUDA toolkit found (v5.0+), disabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
            set(CUDA_NOTFOUND_AUTO ON)
        else()
            # the user requested CUDA, but it wasn't found
            message(FATAL_ERROR "${CUDA_NOTFOUND_MESSAGE}")
        endif()
    else()
        if (GMX_GPU_AUTO)
            message(STATUS "Enabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE ON)
        endif()
    endif() # NOT CUDA_FOUND
endif()

# Try to find NVML if a GPU accelerated binary should be build.
if (GMX_GPU)
    if (DEFINED NVML_LIBRARY)
        set(NVML_FIND_QUIETLY TRUE)
    endif()
    find_package(NVML)
    # TODO Default to off, since linking is not implemented reliably
    option(GMX_USE_NVML "Use NVML support for better CUDA performance" OFF)
    mark_as_advanced(GMX_USE_NVML)
    if(GMX_USE_NVML)
        if(NVML_FOUND)
            include_directories(SYSTEM ${NVML_INCLUDE_DIR})
            set(HAVE_NVML 1)
            list(APPEND GMX_EXTRA_LIBRARIES ${NVML_LIBRARY})
        else()
            message(FATAL_ERROR "NVML support was required, but was not detected. Please consult the install guide.")
        endif()
    endif()
endif()

# Annoyingly enough, FindCUDA leaves a few variables behind as non-advanced.
# We need to mark these advanced outside the conditional, otherwise, if the
# user turns GMX_GPU=OFF after a failed cmake pass, these variables will be
# left behind in the cache.
mark_as_advanced(CUDA_SDK_ROOT_DIR
                 CUDA_USE_STATIC_CUDA_RUNTIME
                 CUDA_dl_LIBRARY CUDA_rt_LIBRARY
                 )
if(NOT GMX_GPU)
    mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
    mark_as_advanced(CUDA_HOST_COMPILER)
endif()

# Try to execute ${CUDA_NVCC_EXECUTABLE} --version and set the output
# (or an error string) in the argument variable.
# Note that semicolon is used as separator for nvcc.
#
# Parameters:
#   COMPILER_INFO   - [output variable] string with compiler path, ID and
#                     some compiler-provided information
#   COMPILER_FLAGS  - [output variable] flags for the compiler
#
macro(get_cuda_compiler_info COMPILER_INFO COMPILER_FLAGS)
    if(NOT GMX_CLANG_CUDA)
        if(CUDA_NVCC_EXECUTABLE)

            # Get the nvcc version string. This is multi-line, but since it is only 4 lines
            # and might change in the future it is better to store than trying to parse out
            # the version from the current format.
            execute_process(COMMAND ${CUDA_NVCC_EXECUTABLE} --version
                RESULT_VARIABLE _nvcc_version_res
                OUTPUT_VARIABLE _nvcc_version_out
                ERROR_VARIABLE  _nvcc_version_err
                OUTPUT_STRIP_TRAILING_WHITESPACE)
            if (${_nvcc_version_res} EQUAL 0)
                # Fix multi-line mess: Replace newline with ";" so we can use it in a define
                string(REPLACE "\n" ";" _nvcc_info_singleline ${_nvcc_version_out})
                SET(${COMPILER_INFO} "${CUDA_NVCC_EXECUTABLE} ${_nvcc_info_singleline}")
                string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
                SET(_compiler_flags "${CUDA_NVCC_FLAGS_${_build_type}}")
                if(CUDA_PROPAGATE_HOST_FLAGS)
                    string(REGEX REPLACE "[ ]+" ";" _cxx_flags_nospace "${BUILD_CXXFLAGS}")
                endif()
                SET(${COMPILER_FLAGS} "${CUDA_NVCC_FLAGS}${CUDA_NVCC_FLAGS_${_build_type}}; ${_cxx_flags_nospace}")
            else()
                SET(${COMPILER_INFO} "N/A")
                SET(${COMPILER_FLAGS} "N/A")
            endif()
        endif()
    else()
        # CXX compiler is the CUDA compiler
        set(${COMPILER_INFO} "${CMAKE_CXX_COMPILER}  ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
        # there are some extra flags
        set(${COMPILER_FLAGS} "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${GMX_CUDA_CLANG_FLAGS}")
    endif()
endmacro ()

macro(enable_multiple_cuda_compilation_units)
    message(STATUS "Enabling multiple compilation units for the CUDA non-bonded module.")
    set_property(CACHE GMX_CUDA_NB_SINGLE_COMPILATION_UNIT PROPERTY VALUE OFF)
endmacro()

include(CMakeDependentOption)
include(gmxOptionUtilities)
macro(gmx_gpu_setup)
    if(GMX_GPU)
        if(NOT GMX_CLANG_CUDA)
            if(NOT CUDA_NVCC_EXECUTABLE)
                message(FATAL_ERROR "nvcc is required for a CUDA build, please set CUDA_TOOLKIT_ROOT_DIR appropriately")
            endif()
            # set up nvcc options
            include(gmxManageNvccConfig)
        else()
            include(gmxManageClangCudaConfig)
        endif()

        gmx_check_if_changed(_cuda_version_changed CUDA_VERSION)

        # Generate CUDA RT API version string which will end up in config.h
        # We do this because nvcc is silly enough to not define its own version
        # (which should match the CUDA runtime API version AFAICT) and we want to
        # avoid creating the fragile dependency on cuda_runtime_api.h.
        #
        # NOTE: CUDA v7.5 is expected to have nvcc define it own version, so in the
        # future we should switch to using that version string instead of our own.
        if (NOT GMX_CUDA_VERSION OR _cuda_version_changed)
            MATH(EXPR GMX_CUDA_VERSION "${CUDA_VERSION_MAJOR}*1000 + ${CUDA_VERSION_MINOR}*10")
        endif()

        if (_cuda_version_changed)
            # check the generated CUDA API version against the one present in cuda_runtime_api.h
            try_compile(_get_cuda_version_compile_res
                ${CMAKE_BINARY_DIR}
                ${CMAKE_SOURCE_DIR}/cmake/TestCUDAVersion.c
                COMPILE_DEFINITIONS "-DGMX_CUDA_VERSION=${GMX_CUDA_VERSION}"
                CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CUDA_TOOLKIT_INCLUDE}"
                OUTPUT_VARIABLE _get_cuda_version_compile_out)

            if (NOT _get_cuda_version_compile_res)
                if (_get_cuda_version_compile_out MATCHES "CUDA version mismatch")
                    message(FATAL_ERROR "The CUDA API version generated internally from the compiler version does not match the version reported by cuda.h. This means either that the CUDA detection picked up mismatching nvcc and the CUDA headers (likely not part of the same toolkit installation) or that there is an error in the internal version generation. If you are sure that it is not the former causing the error (check the relevant cache variables), define the GMX_CUDA_VERSION cache variable to work around the error.")
                else()
                    message(FATAL_ERROR "Could not detect CUDA runtime API version")
                endif()
            endif()
        endif()

        # no OpenMP is no good!
        if(NOT GMX_OPENMP)
            message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading. Without OpenMP a single CPU core can be used with a GPU which is not optimal. Note that with MPI multiple processes can be forced to use a single GPU, but this is typically inefficient. You need to set both C and C++ compilers that support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
        endif()
    endif() # GMX_GPU

    if (GMX_CLANG_CUDA)
        set (_GMX_CUDA_NB_SINGLE_COMPILATION_UNIT_DEFAULT FALSE)
    else()
        set (_GMX_CUDA_NB_SINGLE_COMPILATION_UNIT_DEFAULT TRUE)
    endif()
    cmake_dependent_option(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT
        "Whether to compile the CUDA non-bonded module using a single compilation unit." ${_GMX_CUDA_NB_SINGLE_COMPILATION_UNIT_DEFAULT}
        "GMX_GPU" ON)
    mark_as_advanced(GMX_CUDA_NB_SINGLE_COMPILATION_UNIT)

    if (GMX_GPU AND NOT GMX_CLANG_CUDA)
        # We need to use single compilation unit for kernels:
        # when compiling with nvcc for CC 2.x devices where buggy kernel code is generated
        gmx_check_if_changed(_gmx_cuda_target_changed GMX_CUDA_TARGET_SM GMX_CUDA_TARGET_COMPUTE CUDA_NVCC_FLAGS)

        if(_gmx_cuda_target_changed OR NOT GMX_GPU_DETECTION_DONE)
            if((NOT GMX_CUDA_TARGET_SM AND NOT GMX_CUDA_TARGET_COMPUTE) OR
                (GMX_CUDA_TARGET_SM MATCHES "2[01]" OR GMX_CUDA_TARGET_COMPUTE MATCHES "2[01]"))
                message(STATUS "Enabling single compilation unit for the CUDA non-bonded module. Multiple compilation units are not compatible with CC 2.x devices, to enable the feature specify only CC >=3.0 target architectures in GMX_CUDA_TARGET_SM/GMX_CUDA_TARGET_COMPUTE.")
                set_property(CACHE GMX_CUDA_NB_SINGLE_COMPILATION_UNIT PROPERTY VALUE ON)
            else()
                message(STATUS "Enabling multiple compilation units for the CUDA non-bonded module.")
                set_property(CACHE GMX_CUDA_NB_SINGLE_COMPILATION_UNIT PROPERTY VALUE OFF)
            endif()
        endif()
    endif()
endmacro()
