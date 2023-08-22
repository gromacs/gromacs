#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2023- The GROMACS Authors
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

# This CMake find package follows conventions, namely it sets Hip_FOUND
# cache variable upon success.
# The package creates the shared library target Hip::amdhip
# The package supports following components:
# - hiprtc: create the shared library target Hip::hiprtc

include(FindPackageHandleStandardArgs)

find_library(Hip_amdhip_LIBRARY
    NAMES amdhip64
    PATHS
        "$ENV{HIP_PATH}"
        "$ENV{ROCM_PATH}/hip"
        /opt/rocm/
    PATH_SUFFIXES
        lib
        lib64
    )

if("hiprtc" IN_LIST Hip_FIND_COMPONENTS)
    find_library(Hip_hiprtc_LIBRARY
        NAMES hiprtc
        PATHS
            "$ENV{HIP_PATH}"
            "$ENV{ROCM_PATH}/hip"
            /opt/rocm/
        PATH_SUFFIXES
            lib
            lib64
        )
    find_path(Hip_hiprtc_INCLUDE_DIR
        NAMES hip/hiprtc.h
        PATHS
            "$ENV{HIP_PATH}"
            "$ENV{ROCM_PATH}/hip"
            /opt/rocm/
        PATH_SUFFIXES
            include
            ../include
        )
    if(Hip_hiprtc_LIBRARY AND Hip_hiprtc_INCLUDE_DIR)
        set(Hip_hiprtc_FOUND TRUE)
    endif()
endif()


find_package_handle_standard_args(Hip HANDLE_COMPONENTS REQUIRED_VARS Hip_amdhip_LIBRARY)

if (Hip_FOUND)
    mark_as_advanced(Hip_amdhip_LIBRARY)
    if(NOT TARGET Hip::amdhip)
        add_library(Hip::amdhip SHARED IMPORTED)
        set_property(TARGET Hip::amdhip PROPERTY IMPORTED_LOCATION ${Hip_amdhip_LIBRARY})
        target_compile_definitions(Hip::amdhip INTERFACE -D__HIP_PLATFORM_AMD__)
    endif()
endif()

if(Hip_hiprtc_FOUND)
    mark_as_advanced(Hip_hiprtc_LIBRARY)
    mark_as_advanced(Hip_hiprtc_INCLUDE_DIR)

    if(NOT TARGET Hip::hiprtc)
        add_library(Hip::hiprtc SHARED IMPORTED)
        set_property(TARGET Hip::hiprtc PROPERTY IMPORTED_LOCATION ${Hip_hiprtc_LIBRARY})
        target_include_directories(Hip::hiprtc INTERFACE ${Hip_hiprtc_INCLUDE_DIR})
    endif()
endif()
