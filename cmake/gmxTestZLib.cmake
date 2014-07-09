#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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

# - Define macro to check if linking to zlib actually works,
# because the find_package macro is content if one exists.  This can
# fail in cross-compilation environments, and we want to know about
# zlib so the zlib TNG support is built only when it will work.
#
#  GMX_TEST_ZLIB(VARIABLE)
#
#  VARIABLE will be set to true if zlib support is present

include(CheckLibraryExists)
include(gmxOptionUtilities)
function(GMX_TEST_ZLIB VARIABLE)
    if(NOT ZLIB_FOUND)
        set(${VARIABLE} OFF PARENT_SCOPE)
        return()
    endif()

    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if(${_cmake_build_type} STREQUAL "MSAN")
        # Linking MSan-enabled zlib in a way that can be tested with
        # try_compile setup is tricky, but someone doing an MSan build
        # can take care of themselves.
        set(${VARIABLE} ON PARENT_SCOPE)
        return()
    endif()

    gmx_check_if_changed(_do_zlib_recompile ZLIB_INCLUDE_DIR ZLIB_LIBRARIES)
    if(_do_zlib_recompile)
        unset(ZLIB_LINKS_OK CACHE)
    endif()
    check_library_exists("${ZLIB_LIBRARIES}" "zlibVersion" "" ZLIB_LINKS_OK)
    set(${VARIABLE} ${ZLIB_LINKS_OK} PARENT_SCOPE)
endfunction()



