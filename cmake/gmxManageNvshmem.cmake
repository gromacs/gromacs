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
    
    set(_cuda_arch_message "${GMX_CUDA_ARCHITECTURES}")
    string(SHA1 _cuda_arch_message_hash "${_cuda_arch_message}")
    # Don't report the list of architectures if it's the same as reported before here or in the main CUDA module
    if(NOT LAST_REPORTED_CUDA_ARCH_MESSAGE_HASH_NVSHMEM STREQUAL _cuda_arch_message_hash AND NOT LAST_REPORTED_CUDA_ARCH_MESSAGE_HASH_CUDA STREQUAL _cuda_arch_message_hash)
        message(STATUS "Updated list of CUDA architectures: ${_cuda_arch_message}")
        set(LAST_REPORTED_CUDA_ARCH_MESSAGE_HASH_NVSHMEM "${_cuda_arch_message_hash}" CACHE INTERNAL "The hash of the last reported CUDA architecture list, after NVSHMEM detection")
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(LANGUAGES CUDA RESULT cuda_ipo_supported)
    if (cuda_ipo_supported AND NOT CHECK_CUDA_IPO_QUIETLY)
        message(STATUS "CUDA IPO Supported: ${cuda_ipo_supported}")
        set(CHECK_CUDA_IPO_QUIETLY TRUE CACHE INTERNAL "Be quiet about CUDA IPO support")
    endif()

    # Helper function to create NVSHMEM separable compilation libraries
    # Why separate static libraries?
    # - Main reason: avoid enabling CUDA_SEPARABLE_COMPILATION globally. Global
    #   separable compilation can lead to suboptimal kernel performance compared
    #   to whole-program compilation. By restricting separable compilation to
    #   only the NVSHMEM-using modules (per-module static libs), we preserve
    #   optimal codegen for the rest of the CUDA sources.
    # - NVSHMEM provides both host and device components. Device code that calls
    #   NVSHMEM still needs separable compilation so that all device objects are
    #   collected and device-linked correctly.
    # - Building per-module static libraries with CUDA_SEPARABLE_COMPILATION=ON
    #   lets us aggregate device code across directories and only run device link
    #   where required.
    # - Some toolchains benefit from hiding device symbols during device link to
    #   avoid multiple-definition issues (see -fvisibility=hidden used with
    #   $<DEVICE_LINK:...> below).
    # When to use:
    # - Whenever adding CUDA device code that calls NVSHMEM APIs.
    function(gmx_add_library_for_nvshmem)
        cmake_parse_arguments(ARG
            "STATIC;SHARED"
            "TARGET_NAME"
            "SOURCES;PUBLIC_TARGET_DEPENDENCIES;PRIVATE_TARGET_DEPENDENCIES;INTERFACE_TARGET_DEPENDENCIES"
            ${ARGN}
        )

        if(NOT ARG_TARGET_NAME)
            message(FATAL_ERROR "gmx_add_library_for_nvshmem: TARGET_NAME is required")
        endif()
        if(NOT ARG_SOURCES)
            message(FATAL_ERROR "gmx_add_library_for_nvshmem: SOURCES is required")
        endif()
        if(ARG_STATIC AND ARG_SHARED)
            message(FATAL_ERROR "gmx_add_library_for_nvshmem: Cannot have both STATIC and SHARED")
        endif()
        # Default to STATIC if no library type specified
        set(LIB_TYPE STATIC)
        if(ARG_SHARED)
            set(LIB_TYPE SHARED)
        endif()

        # Create the separable compilation library
        add_library(${ARG_TARGET_NAME} ${LIB_TYPE} ${ARG_SOURCES})
        # Set CUDA separable compilation properties
        set_target_properties(${ARG_TARGET_NAME} PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_RESOLVE_DEVICE_SYMBOLS OFF
            POSITION_INDEPENDENT_CODE ON
            CUDA_ARCHITECTURES "${GMX_CUDA_ARCHITECTURES}")
        # Apply standard GROMACS compile options
        gmx_target_compile_options(${ARG_TARGET_NAME})
        target_compile_definitions(${ARG_TARGET_NAME} PRIVATE HAVE_CONFIG_H TMPI_EXPORTS PUBLIC TMPI_USE_VISIBILITY)
        target_compile_options(${ARG_TARGET_NAME} PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:${GMX_CUDA_FLAGS}>")
        # NVSHMEM-specific linking options
        # -fvisibility=hidden must be passed during the device link of ${ARG_TARGET_NAME}
        # to hide libnvshmem_device.a symbols from the ${ARG_TARGET_NAME} device link output
        # that is later aggregated into libgromacs. This allows applications/tests to also link
        # libnvshmem_device.a for their own NVSHMEM kernels while linking against libgromacs.
        # Omitting this causes runtime failures due to conflicting device symbols.
        target_link_options(${ARG_TARGET_NAME} PRIVATE $<DEVICE_LINK:-fvisibility=hidden>)
        target_link_libraries(${ARG_TARGET_NAME} PRIVATE nvshmem_device_lib nvshmem_host_lib)
        # Standard GROMACS include directories
        target_include_directories(${ARG_TARGET_NAME} PRIVATE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src/include>)
        # Enable device LTO if supported
        if(cuda_ipo_supported)
            set_target_properties(${ARG_TARGET_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
        endif()
        # Link target dependencies
        if(ARG_PUBLIC_TARGET_DEPENDENCIES)
            target_link_libraries(${ARG_TARGET_NAME} PUBLIC ${ARG_PUBLIC_TARGET_DEPENDENCIES})
        endif()
        if(ARG_PRIVATE_TARGET_DEPENDENCIES)
            target_link_libraries(${ARG_TARGET_NAME} PRIVATE ${ARG_PRIVATE_TARGET_DEPENDENCIES})
        endif()
        if(ARG_INTERFACE_TARGET_DEPENDENCIES)
            target_link_libraries(${ARG_TARGET_NAME} INTERFACE ${ARG_INTERFACE_TARGET_DEPENDENCIES})
        endif()
    endfunction()
endif()
