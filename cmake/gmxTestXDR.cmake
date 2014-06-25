#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2014, by the GROMACS development team, led by
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

# - Define macro to check if system XDR routines are present
#
#  GMX_TEST_XDR(VARIABLE)
#
#  VARIABLE will be set to true if XDR support is present
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_XDR VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for system XDR support")

	# First check without any special flags
        TRY_COMPILE(XDR_COMPILE_OK "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestXDR.c")

        if(XDR_COMPILE_OK)
	    MESSAGE(STATUS "Checking for system XDR support - present")
        else()
            MESSAGE(STATUS "Checking for system XDR support - not present")
      	endif()

        set(${VARIABLE} ${XDR_COMPILE_OK} CACHE INTERNAL "Result of test for system XDR support" FORCE)

    ENDIF()
ENDMACRO(GMX_TEST_XDR VARIABLE)



