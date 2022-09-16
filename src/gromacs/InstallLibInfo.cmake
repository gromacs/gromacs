#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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

function (do_pkgconfig)
    set(PKG_CFLAGS "")
    foreach (_dir ${INSTALLED_HEADER_INCLUDE_DIRS})
        if (IS_ABSOLUTE ${_dir})
            set(PKG_CFLAGS "${PKG_CFLAGS} -I${_dir}")
        else()
            set(PKG_CFLAGS "${PKG_CFLAGS} -I${CMAKE_INSTALL_PREFIX}/${_dir}")
        endif()
    endforeach()
    if (INSTALLED_HEADER_DEFINITIONS)
        foreach (_def ${INSTALLED_HEADER_DEFINITIONS})
            set(PKG_CFLAGS "${PKG_CFLAGS} ${_def}")
        endforeach()
    endif()
    set(PKG_CFLAGS "${PKG_CFLAGS} ${OpenMP_C_FLAGS}")

    configure_file(libgromacs.pc.cmakein
                   libgromacs.pc @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libgromacs.pc
            DESTINATION ${GMX_INSTALL_PKGCONFIGDIR}
            RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
            COMPONENT development)
endfunction()

function (do_cmake_config)
    # Install everything into a subdirectory, because
    #  1. CMake expects things to be there for CMAKE_PREFIX_PATH to work, and
    #  2. This nicely isolates files for different suffixes from each other.
    set(GMX_INSTALL_CMAKEPKGDIR ${GMX_INSTALL_CMAKEDIR}/gromacs${GMX_LIBS_SUFFIX})

    # Install import definitions that take care of the library locations and
    # library dependencies.
    set(EXPORT_FILE_NAME libgromacs.cmake)
    if (NOT BUILD_SHARED_LIBS)
        set(EXPORT_FILE_NAME libgromacs_static.cmake)
    endif()
    install(EXPORT libgromacs
            FILE ${EXPORT_FILE_NAME}
            NAMESPACE Gromacs::
            DESTINATION ${GMX_INSTALL_CMAKEPKGDIR}
            COMPONENT libraries)

    get_filename_component(GROMACS_CXX_COMPILER ${CMAKE_CXX_COMPILER} REALPATH)
    if (CMAKE_OSX_DEPLOYMENT_TARGET OR CMAKE_OSX_ARCHITECTURES)
        set(_gmx_osx_config
"SET(CMAKE_OSX_DEPLOYMENT_TARGET \"${CMAKE_OSX_DEPLOYMENT_TARGET}\" CACHE STRING \"GROMACS Deployment target.\"
SET(CMAKE_OSX_ARCHITECTURES \"${CMAKE_OSX_ARCHITECTURES}\" CACHE STRING \"GROMACS architectures.\")")
    endif ()

    if (GMX_LIB_MPI)
        set(_gmx_mpi_config
"SET(MPI_C_COMPILER \"${MPI_C_COMPILER}\")
SET(MPI_CXX_COMPILER \"${MPI_CXX_COMPILER}\")")
    endif ()
    configure_file(gromacs-config.cmake.cmakein
                   gromacs-config.cmake @ONLY)
    configure_file(gromacs-config-version.cmake.cmakein
                   gromacs-config-version.cmake @ONLY)
    configure_file(gromacs-hints.in.cmake
                   gromacs-hints.cmake @ONLY)
    unset(_gmx_mpi_config)
    unset(_gmx_osx_config)
    option(GMX_REQUIRE_VALID_CMAKE_HINTS "Force CMake error if generated hints are not usable." OFF)
    mark_as_advanced(GMX_REQUIRE_VALID_CMAKE_HINTS)
    if (GMX_REQUIRE_VALID_CMAKE_HINTS)
        # Test the generated toolchain file.
        set(TEMPDIR "${CMAKE_CURRENT_BINARY_DIR}/cmake-configure-test")
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake-configure-test)
        execute_process(COMMAND
                            ${CMAKE_COMMAND}
                            -G "${CMAKE_GENERATOR}"
                            -C ${CMAKE_CURRENT_BINARY_DIR}/gromacs-hints.cmake
                            -DGMX_REQUIRE_VALID_CMAKE_HINTS=FALSE
                            ${CMAKE_SOURCE_DIR}
                        RESULT_VARIABLE result
                        OUTPUT_VARIABLE output
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        ERROR_VARIABLE output
                        ERROR_STRIP_TRAILING_WHITESPACE
                        WORKING_DIRECTORY ${TEMPDIR})
        if (result)
            message(FATAL_ERROR "Generated gromacs-hints.cmake does not produce a valid CMake environment: ${output}")
        else()
            # Note: We confirm that the hints file does not result in errors, but we do not check
            # the degree to which it produces a compatible build system. As inadequacies are
            # encountered, we can add additional checks here.
            message(STATUS "Verified usability of gromacs-hints.cmake")
            # We clean up after ourselves
            FILE(REMOVE_RECURSE ${TEMPDIR})
        endif ()
    endif ()

    # The configuration files are also installed with the suffix, even though
    # the directory already contains the suffix. This allows simple
    # find_package(GROMACS NAMES gromacs_d) to find them, without also
    # specifying CONFIGS.
    # Note, however, that the exports file for Gromacs::libgromacs exists with
    # the same name in each of the gromacs${GMX_LIBS_SUFFIX} subdirectories.
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gromacs-config.cmake
            DESTINATION ${GMX_INSTALL_CMAKEPKGDIR}
            RENAME "gromacs${GMX_LIBS_SUFFIX}-config.cmake"
            COMPONENT development)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gromacs-config-version.cmake
            DESTINATION ${GMX_INSTALL_CMAKEPKGDIR}
            RENAME "gromacs${GMX_LIBS_SUFFIX}-config-version.cmake"
            COMPONENT development)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gromacs-hints.cmake
            DESTINATION ${GMX_INSTALL_CMAKEPKGDIR}
            RENAME "gromacs-hints${GMX_LIBS_SUFFIX}.cmake"
            COMPONENT development)
endfunction()

do_pkgconfig()
do_cmake_config()
