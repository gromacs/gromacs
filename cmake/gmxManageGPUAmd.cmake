#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
# - search for ROCM and set GMX_GPU=ON we find it
# - check whether GPUs are present
# - if ROCM is not found but GPUs were detected issue a warning

# detect GPUs in the build host machine
if ((GMX_GPU OR GMX_GPU_AUTO) AND NOT GMX_GPU_DETECTION_DONE)
    include(gmxDetectGpuAmd)
    gmx_detect_gpu()
endif()

# We need to call find_package even when we've already done the detection/setup
if( (GMX_GPU OR GMX_GPU_AUTO) AND GMX_DETECT_GPU_AVAILABLE)
    find_package(hcc QUIET CONFIG PATHS /opt/rocm)
    find_package(hip QUIET CONFIG PATHS /opt/rocm)

    if ( hip_FOUND ) 
         message(STATUS "HIP found in gmxManageGPU !")
    endif()
endif()

if((GMX_GPU OR GMX_GPU_AUTO) AND NOT GMX_GPU_DETECTION_DONE)
    # assemble warning/error message
    if (GMX_DETECT_GPU_AVAILABLE)
        set(_msg "${GMX_DETECT_GPU_COUNT} ROCM GPU(s) found in the system")

        # append GPU names
        if (NOT GMX_DETECT_GPU_INFO STREQUAL "")
            set(_msg "${_msg}:")
            foreach(gpu ${GMX_DETECT_GPU_INFO})
                set(_msg "${_msg} ${gpu} ")
            endforeach()
        endif()
    endif()

    set(ROCM_NOTFOUND_MESSAGE "mdrun supports native GPU acceleration on ROCM hardward). This requires the ROCM HIP API, which was not found. The typical location would be /opt/rocm. Note that CPU or GPU acceleration can be selected at runtime.  ${_msg}")
    unset(_msg)

    if (NOT hip_FOUND)
        if (GMX_GPU_AUTO)
            # Disable GPU acceleration in auto mode
            message(STATUS "No compatible ROCM found, disabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE OFF)
            set(ROCM_NOTFOUND_AUTO ON)
        else()
            # the user requested ROCM, but it wasn't found
            message(FATAL_ERROR "${ROCM_NOTFOUND_MESSAGE}")
        endif()
    else()
        if (GMX_GPU_AUTO)
            message(STATUS "Enabling native GPU acceleration")
            set_property(CACHE GMX_GPU PROPERTY VALUE ON)
        endif()
    endif()
endif()

macro(get_hip_compiler_info COMPILER_INFO COMPILER_FLAGS)
    find_program(HIP_CONFIG hipconfig
         PATH_SUFFIXES bin
         PATHS /opt/rocm/hip
    )
    if(HIP_CONFIG)
        execute_process(COMMAND ${HIP_CONFIG} --version
                RESULT_VARIABLE _hipcc_version_res
                OUTPUT_VARIABLE _hipcc_version_out
                ERROR_VARIABLE  _hipcc_version_err
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if (${_hipcc_version_res} EQUAL 0)
            find_program(HIP_CC_COMPILER hipcc
                 PATH_SUFFIXES bin
                 PATHS /opt/rocm/hip
            )
            if (HIP_CC_COMPILER)
                SET(${COMPILER_INFO} "${HIP_CC_COMPILER} ${_hipcc_version_out}")
                SET(${COMPILER_FLAGS} "")
            else()
                SET(${COMPILER_INFO} "N/A")
                SET(${COMPILER_FLAGS} "N/A")
            endif()
        else()
            SET(${COMPILER_INFO} "N/A")
            SET(${COMPILER_FLAGS} "N/A")
        endif()
        # there are some extra flags
        set(${COMPILER_FLAGS} "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}}")
     endif()
endmacro ()

macro(enable_multiple_rocm_compilation_units)
    message(STATUS "Enabling multiple compilation units for the ROCM non-bonded module.")
    set_property(CACHE GMX_ROCM_NB_SINGLE_COMPILATION_UNIT PROPERTY VALUE OFF)
endmacro()

include(CMakeDependentOption)
include(gmxOptionUtilities)
macro(gmx_gpu_setup)
    if(GMX_GPU) 
       find_package(rocfft QUIET CONFIG PATHS /opt/rocm )
       if (NOT rocfft_FOUND ) 
           message(FATAL_ERROR "rocfft is required, but it is not found on this building environment") 
       else() 
           message(STATUS "rocfft is found!") 
       endif() 
    endif() 

    if(GMX_GPU)
        # no OpenMP is no good!
        if(NOT GMX_OPENMP)
            message(WARNING "To use GPU acceleration efficiently, mdrun requires OpenMP multi-threading. Without OpenMP a single CPU core can be used with a GPU which is not optimal. Note that with MPI multiple processes can be forced to use a single GPU, but this is typically inefficient. You need to set both C and C++ compilers that support OpenMP (CC and CXX environment variables, respectively) when using GPUs.")
        endif()
    endif() # GMX_GPU

    option(GMX_ROCM_NB_SINGLE_COMPILATION_UNIT "Whether to compile the ROCM non-bonded module using a single compilation unit." OFF)
    mark_as_advanced(GMX_ROCM_NB_SINGLE_COMPILATION_UNIT)

endmacro()
