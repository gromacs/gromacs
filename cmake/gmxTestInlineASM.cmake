#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2009- The GROMACS Authors
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

# - Define macro to check GCC x86 inline ASM support
#
#  GMX_TEST_INLINE_ASM_GCC_X86(VARIABLE)
#
#  VARIABLE will be set to true if GCC x86 inline asm works.

MACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for GCC x86 inline asm")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_gcc_x86.cpp"
                    OUTPUT_VARIABLE INLINE_ASM_COMPILE_OUTPUT)

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for GCC x86 inline asm - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
        else()
            MESSAGE(STATUS "Checking for GCC x86 inline asm - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
        endif()

    ENDIF()
ENDMACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)




