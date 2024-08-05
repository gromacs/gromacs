#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2024- The GROMACS Authors
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

if(GMX_DOUBLE)
    message(FATAL_ERROR "HIP acceleration is not available in double precision")
endif()

set(CMAKE_HIP_STANDARD ${CMAKE_CXX_STANDARD})
set(CMAKE_HIP_STANDARD_REQUIRED ON)

# Using the required version directly doesn't work due to the way the versioning is implemented in HIP
find_package(HIP REQUIRED CONFIG PATHS $ENV{ROCM_PATH} "/opt/rocm")
if (${HIP_VERSION} VERSION_LESS ${REQUIRED_HIP_VERSION})
    message(FATAL_ERROR "The found HIP version ${HIP_VERSION} is less than the required version ${REQUIRED_HIP_VERSION}. Please update your ROCm stack")
endif()

enable_language(HIP)
set(GMX_GPU_HIP ON)


find_package(rocprim REQUIRED CONFIG HINTS ${HIP_PACKAGE_PREFIX_DIR})

if(GMX_GPU_FFT_VKFFT)
    include(gmxManageVkFft)
elseif(GMX_GPU_FFT_HIPFFT OR GMX_GPU_FFT_ROCFFT OR GMX_USE_Heffte)
    if (GMX_GPU_FFT_HIPFFT)
        find_package(hipfft REQUIRED CONFIG HINTS ${HIP_PACKAGE_PREFIX_DIR})
    endif()
    find_package(rocfft REQUIRED CONFIG HINTS ${HIP_PACKAGE_PREFIX_DIR})
else()
    message(FATAL_ERROR "The configured GPU FFT library ${GMX_GPU_FFT_LIBRARY} can not be used together with the HIP backend") 
endif()

macro(get_hip_compiler_info COMPILER_INFO COMPILER_FLAGS)
    if(HIP_HIPCONFIG_EXECUTABLE)
        execute_process(COMMAND ${HIP_HIPCONFIG_EXECUTABLE} --version
                RESULT_VARIABLE _hipcc_version_res
                OUTPUT_VARIABLE _hipcc_version_out
                ERROR_VARIABLE  _hipcc_version_err
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(${_hipcc_version_res} EQUAL 0)
            if(HIP_HIPCC_EXECUTABLE)
                set(${COMPILER_INFO} "${HIP_HIPCC_EXECUTABLE} ${_hipcc_version_out}")
                set(${COMPILER_FLAGS} "${GMX_HIP_HIPCC_FLAGS}")
            else()
                message(WARNING "The HIP_HIPCC_EXECUTABLE variable has not been set, so no compiler information is available")
                set(${COMPILER_INFO} "N/A")
                set(${COMPILER_FLAGS} "N/A")
            endif()
        else()
            set(${COMPILER_INFO} "N/A")
            set(${COMPILER_FLAGS} "N/A")
        endif()
    else()
        message(WARNING "The HIP_HIPCONFIG_EXECUTABLE variable has not been set by find_package(HIP), meaning GROMACS can not properly present the information about the ROCm toolkit")
    endif()
endmacro()

set(GMX_HIPCC_EXTRA_FLAGS "" CACHE STRING "Extra GROMACS specific HIPCC compiler flags")
include(gmxManageHipccConfig)

option(GMX_HIP_NB_SINGLE_COMPILATION_UNIT "Whether to compile the HIP non-bonded module using a single compilation unit." OFF)
mark_as_advanced(GMX_HIP_NB_SINGLE_COMPILATION_UNIT)
