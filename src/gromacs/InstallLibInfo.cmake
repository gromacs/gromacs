#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014, by the GROMACS development team, led by
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
            DESTINATION ${PKGCONFIG_INSTALL_DIR}
            RENAME "libgromacs${GMX_LIBS_SUFFIX}.pc"
            COMPONENT development)
endfunction()

function (do_cmake_config)
    # Install everything into a subdirectory, because
    #  1. CMake expects things to be there for CMAKE_PREFIX_PATH to work, and
    #  2. This nicely isolates files for different suffixes from each other.
    set(CMAKE_PACKAGE_DIR ${CMAKE_INSTALL_DIR}/gromacs${GMX_LIBS_SUFFIX})

    # Install import definitions that take care of the library locations and
    # library dependencies.
    set(EXPORT_FILE_NAME libgromacs.cmake)
    if (NOT BUILD_SHARED_LIBS)
        set(EXPORT_FILE_NAME libgromacs_static.cmake)
    endif()
    install(EXPORT libgromacs
            FILE ${EXPORT_FILE_NAME}
            DESTINATION ${CMAKE_PACKAGE_DIR}
            COMPONENT libraries)

    get_filename_component(GROMACS_CXX_COMPILER ${CMAKE_CXX_COMPILER} REALPATH)
    configure_file(gromacs-config.cmake.cmakein
                   gromacs-config.cmake @ONLY)
    configure_file(gromacs-config-version.cmake.cmakein
                   gromacs-config-version.cmake @ONLY)
    # The configuration files are also installed with the suffix, even though
    # the directory already contains the suffix. This allows simple
    # find_package(GROMACS NAMES gromacs_d) to find them, without also
    # specifying CONFIGS.
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gromacs-config.cmake
            DESTINATION ${CMAKE_PACKAGE_DIR}
            RENAME "gromacs${GMX_LIBS_SUFFIX}-config.cmake"
            COMPONENT development)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gromacs-config-version.cmake
            DESTINATION ${CMAKE_PACKAGE_DIR}
            RENAME "gromacs${GMX_LIBS_SUFFIX}-config-version.cmake"
            COMPONENT development)
endfunction()

do_pkgconfig()
do_cmake_config()
