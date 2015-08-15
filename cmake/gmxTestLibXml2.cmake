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

# - Define macro to check if linking to libxml2 actually works,
# because the find_package macro is content if one exists.  This can
# fail in cross-compilation environments, and we want to know about
# libxml2 so the test binaries are built only when they will work.
#
#  GMX_TEST_LIBXML2(VARIABLE)
#
#  VARIABLE will be set to true if libxml2 support is present

include(CheckLibraryExists)
include(CheckIncludeFiles)
include(gmxOptionUtilities)
function(GMX_TEST_LIBXML2 VARIABLE)
    if(NOT LIBXML2_FOUND)
        set(${VARIABLE} OFF PARENT_SCOPE)
	return()
    endif()

    if(HAVE_ZLIB)
        set(LIBXML2_LIBRARIES "${LIBXML2_LIBRARIES};${ZLIB_LIBRARIES}" PARENT_SCOPE) #not needed for dynamic but does not hurt
    endif()

    string(TOUPPER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    if(${_cmake_build_type} STREQUAL "MSAN")
        # Linking MSan-enabled libxml2 and zlib in a way that can be
        # tested with try_compile setup is tricky, but someone doing
        # an MSan build can take care of themselves (for now, at
        # least).
        set(${VARIABLE} ON PARENT_SCOPE)
        return()
    endif()

    gmx_check_if_changed(_do_libxml2_recompile LIBXML2_INCLUDE_DIR LIBXML2_LIBRARIES)
    if(_do_libxml2_recompile)
        unset(LIBXML2_LINKS_OK CACHE)
    endif()
    check_library_exists("${LIBXML2_LIBRARIES}" "xmlTextWriterEndAttribute" "" LIBXML2_LINKS_OK)
    if(LIBXML2_LINKS_OK)
        #check that xml headers can be included
        set(CMAKE_REQUIRED_INCLUDES "${LIBXML2_INCLUDE_DIR}")
	check_include_files("libxml/parser.h" LIBXML2_INCL_OK)
	if(NOT LIBXML2_INCL_OK)
            #xml headers depend on iconv.h. Test whether adding its path fixes the problem
            find_path(ICONV_INCLUDE_DIR iconv.h)
            if(ICONV_INCLUDE_DIR)
                set(CMAKE_REQUIRED_INCLUDES "${LIBXML2_INCLUDE_DIR};${ICONV_INCLUDE_DIR}")
		unset(LIBXML2_INCL_OK CACHE)
		check_include_files("libxml/parser.h" LIBXML2_INCL_OK)
		set(LIBXML2_INCLUDE_DIR "${LIBXML2_INCLUDE_DIR};${ICONV_INCLUDE_DIR}" CACHE PATH "Libxml2 include path" FORCE)
            endif()
        endif()
	set(${VARIABLE} ${LIBXML2_INCL_OK} PARENT_SCOPE)
    else()
        set(${VARIABLE} OFF PARENT_SCOPE)
    endif()
endfunction()



