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

if(GMX_NVSHMEM)
    if(NOT GMX_GPU_CUDA)
        message(FATAL_ERROR "NVSHMEM support requires a CUDA build")
    endif()
    if(NOT GMX_LIB_MPI)
        message(FATAL_ERROR "NVSHMEM support requires a library MPI build")
    endif()
    if (GMX_CLANG_CUDA)
        message(FATAL_ERROR "NVSHMEM is not supported with Clang CUDA build")
    endif()
    if (GMX_USE_CUFFTMP)
        message(FATAL_ERROR "Direct use of NVSHMEM is not yet supported together with cuFFTMp (which uses NVSHMEM internally). GMX_NVSHMEM and GMX_USE_CUFFTMP cannot be enabled at the same time.")
    endif()

    find_library(NVSHMEM_DEVICE_LIBS NAMES nvshmem_device PATHS "${GMX_NVSHMEM_HOME}/lib/" REQUIRED)
    find_library(NVSHMEM_HOST_LIBS NAMES nvshmem_host PATHS "${GMX_NVSHMEM_HOME}/lib/" REQUIRED)
    find_path(NVSHMEM_INCLUDE NAMES nvshmem.h PATHS "${GMX_NVSHMEM_HOME}/include/" REQUIRED)
    add_library(nvshmem_host_lib SHARED IMPORTED GLOBAL)
    set_target_properties(nvshmem_host_lib PROPERTIES IMPORTED_LOCATION ${NVSHMEM_HOST_LIBS})
    set_target_properties(nvshmem_host_lib PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES CUDA)
    target_include_directories(nvshmem_host_lib INTERFACE $<BUILD_INTERFACE:${NVSHMEM_INCLUDE}>)
    target_link_libraries(nvshmem_host_lib INTERFACE CUDA::nvml CUDA::cuda_driver)

    add_library(nvshmem_device_lib STATIC IMPORTED GLOBAL)
    set_target_properties(nvshmem_device_lib PROPERTIES IMPORTED_LOCATION ${NVSHMEM_DEVICE_LIBS})
    target_include_directories(nvshmem_device_lib INTERFACE $<BUILD_INTERFACE:${NVSHMEM_INCLUDE}>)
    set_target_properties(nvshmem_device_lib PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES CUDA)

    # Since NVSHMEM 3.06, minium device support is Volta (SM 70+) so
    # we filter all the archs below SM 70 from GMX_CUDA_ARCHITECTURES
    foreach(_old_arch 35 37 50 52 53 60 61 62)
        list(REMOVE_ITEM GMX_CUDA_ARCHITECTURES "${_old_arch}" "${_old_arch}-real" "${_old_arch}-virtual")
    endforeach()
    message(STATUS "Updated list of CUDA architectures: ${GMX_CUDA_ARCHITECTURES}")
endif()
