#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

# CMake 3.0-3.1 has a bug in the following case, which breaks
# configuration on at least BlueGene/Q. Fixed in 3.1.1
if ((NOT CMAKE_VERSION VERSION_LESS "3.0.0") AND
    (CMAKE_VERSION VERSION_LESS "3.1.1") AND
        (CMAKE_CROSSCOMPILING AND NOT CMAKE_SYSTEM_PROCESSOR))
    message(STATUS "Cannot search for CUDA because the CMake find package has a bug. Set a valid CMAKE_SYSTEM_PROCESSOR if you need to detect CUDA")
else()
    set(CAN_RUN_CUDA_FIND_PACKAGE 1)
endif()

# We need to call find_package even when we've already done the detection/setup
if(GMX_GPU OR GMX_GPU_AUTO AND CAN_RUN_CUDA_FIND_PACKAGE)
    if(NOT GMX_GPU AND NOT GMX_DETECT_GPU_AVAILABLE)
        # Stay quiet when detection has occured and found no GPU.
        # Noise is acceptable when there is a GPU or the user required one.
        set(FIND_CUDA_QUIETLY QUIET)
    endif()
    find_package(CUDA ${REQUIRED_CUDA_VERSION} ${FIND_CUDA_QUIETLY})

    # The IBM xlc compiler chokes if we use both altivec and Cuda. Solve
    # this by not propagating the flags in this case, but add -O3
    # to make sure we don't turn off optimization.
    if(CMAKE_CXX_COMPILER_ID MATCHES "XL")
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
        list(APPEND CUDA_NVCC_FLAGS "-O3")
    endif()

    # Cmake 2.8.12 (and CMake 3.0) introduced a new bug where the cuda
    # library dir is added twice as an rpath on APPLE, which in turn causes
    # the install_name_tool to wreck the binaries when it tries to remove this
    # path. Since this is set inside the cuda module, we remove the extra rpath
    # added in the library string - an rpath is not a library anyway, and at
    # least for Gromacs this works on all CMake versions. This should be
    # reasonably future-proof, since newer versions of CMake appear to handle
    # the rpath automatically based on the provided library path, meaning
    # the explicit rpath specification is no longer needed.
    if(APPLE AND (CMAKE_VERSION VERSION_GREATER 2.8.11))
        foreach(elem ${CUDA_LIBRARIES})
            if(elem MATCHES "-Wl,.*")
                list(REMOVE_ITEM CUDA_LIBRARIES ${elem})
            endif()
        endforeach(elem)
    endif()
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
            message(STATUS "No compatible CUDA toolkit found (v4.0+), disabling native GPU acceleration")
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
    find_package(NVML)
    if(NVML_FOUND)
        include_directories(${NVML_INCLUDE_DIR})
        set(HAVE_NVML 1)
        list(APPEND GMX_EXTRA_LIBRARIES ${NVML_LIBRARY})
    endif(NVML_FOUND)
endif()

# Annoyingly enough, FindCUDA leaves a few variables behind as non-advanced.
# We need to mark these advanced outside the conditional, otherwise, if the
# user turns GMX_GPU=OFF after a failed cmake pass, these variables will be
# left behind in the cache.
mark_as_advanced(CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD)
if(NOT GMX_GPU)
    mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
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
endmacro ()

macro(gmx_gpu_setup)
    # set up nvcc options
    include(gmxManageNvccConfig)

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

    # texture objects are supported in CUDA 5.0 and later
    if (CUDA_VERSION VERSION_GREATER 4.999)
        set(HAVE_CUDA_TEXOBJ_SUPPORT 1)
    endif()

    # Atomic operations used for polling wait for GPU
    # (to avoid the cudaStreamSynchronize + ECC bug).
    # ThreadMPI is now always included. Thus, we don't check for Atomics anymore here.

    # no OpenMP is no good!
    if(NOT GMX_OPENMP)
        message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading. Without OpenMP a single CPU core can be used with a GPU which is not optimal. Note that with MPI multiple processes can be forced to use a single GPU, but this is typically inefficient. You need to set both C and C++ compilers that support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
    endif()
endmacro()
