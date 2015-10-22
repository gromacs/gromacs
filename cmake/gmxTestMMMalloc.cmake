#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2012,2013,2014,2015, by the GROMACS development team, led by
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

#
#  GMX_TEST_MM_MALLOC(VARIABLE)
#
#  VARIABLE will be set to true if we find _mm_malloc() and _mm_free().
#
#  The macro will also check whether the headers mm_malloc.h, malloc.h and
#  xmmintrin.h are present, since the routines might be defined in either,
#  and set the corresponding HAVE_MM_MALLOC_H, HAVE_MALLOC_H, and
#  HAVE_XMMINTRIN_H if it hasn't already been done outside this file.
#
MACRO(GMX_TEST_MM_MALLOC VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        check_include_files(mm_malloc.h HAVE_MM_MALLOC_H)
        check_include_files(malloc.h    HAVE_MALLOC_H)
        check_include_files(xmmintrin.h HAVE_XMMINTRIN_H)

        MESSAGE(STATUS "Checking for _mm_malloc()")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestMMMalloc.cpp"
                    COMPILE_DEFINITIONS
                    "-DHAVE_MM_MALLOC_H=${HAVE_MM_MALLOC_H}"
                    "-DHAVE_MALLOC_H=${HAVE_MALLOC_H}"
                    "-DHAVE_XMMINTRIN_H=${HAVE_XMMINTRIN_H}"
                    OUTPUT_VARIABLE MM_MALLOC_COMPILE_OUTPUT)

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for _mm_malloc() - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for _mm_malloc" FORCE)
        else()
            MESSAGE(STATUS "Checking for _mm_malloc() - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for _mm_malloc()" FORCE)
        endif()

    ENDIF()
ENDMACRO(GMX_TEST_MM_MALLOC VARIABLE)




