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

# Find libxml2 for use in GROMACS. Currently nothing uses it, so this
# file is never included.
#
# The find_package() exports LIBXML2_FOUND, which we should not use
# because it does not tell us that linking will succeed. Instead, we
# test that next.

# - Define macro to check if linking to libxml2 actually works,
# because the find_package macro is content if one exists.  This can
# fail in cross-compilation environments, and we want to know about
# libxml2 so the binaries are built only when they will work.
#
#  GMX_TEST_LIBXML2(VARIABLE)
#
#  VARIABLE will be set to true if libxml2 support is present

include(CheckLibraryExists)
include(gmxOptionUtilities)
function(GMX_TEST_LIBXML2 VARIABLE)
    if(LIBXML2_FOUND)
        gmx_check_if_changed(_do_libxml2_recompile LIBXML2_INCLUDE_DIR LIBXML2_LIBRARIES)
        if(_do_libxml2_recompile)
            unset(LIBXML2_LINKS_OK CACHE)
        endif()
        check_library_exists("${LIBXML2_LIBRARIES}" "xmlTextWriterEndAttribute" "" LIBXML2_LINKS_OK)
        set(${VARIABLE} ${LIBXML2_LINKS_OK} PARENT_SCOPE)
    else()
        set(${VARIABLE} OFF PARENT_SCOPE)
    endif()
endfunction()


# Stay quiet if detection has already run
if(DEFINED LIBXML2_LIBRARIES)
  set(LibXml2_FIND_QUIETLY TRUE)
endif()
find_package(LibXml2)

# Test linking really works
gmx_test_libxml2(HAVE_LIBXML2)

option(GMX_XML "Use libxml2 to parse xml files (currently has no effect)" ${HAVE_LIBXML2})
set(PKG_XML "")
mark_as_advanced(GMX_XML)

if(GMX_XML AND NOT HAVE_LIBXML2)
    message(FATAL_ERROR "libxml2 not found. Set GMX_XML=OFF to compile without XML support")
endif()
if(GMX_XML)
    include_directories(${LIBXML2_INCLUDE_DIR})
    set(PKG_XML libxml-2.0)
    set(XML_LIBRARIES ${LIBXML2_LIBRARIES})
    # If ever using this again, arrange to add dependency on these
    # libraries to the appropriate linking infrastructure
endif()
