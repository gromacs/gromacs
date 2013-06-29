#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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

if(NOT GMX_CPU_ACCELERATION STREQUAL "NONE")
    gmx_add_unit_test(SIMDTests simd-test
        utils.cpp
        load_table.cpp
        load_exclusions.cpp
        arithmetic.cpp
        logic.cpp
        masking.cpp
        load.cpp
        set.cpp
        store.cpp
        )

    # Suppress warnings when the template value parameter for the
    # number of input values has size 0. That value leads to correct
    # behaviour in for loops, i.e. not looping over the number of
    # input values.
    if(CMAKE_COMPILER_IS_GNUCC AND C_COMPILER_VERSION AND C_COMPILER_VERSION VERSION_LESS 4.7)
        set_property(SOURCE set.cpp PROPERTY COMPILE_FLAGS -Wno-type-limits)
    endif()
    if (${CMAKE_C_COMPILER_ID} MATCHES "Clang")
        set_property(SOURCE set.cpp PROPERTY COMPILE_FLAGS -Wno-tautological-compare)
    endif()
endif()
