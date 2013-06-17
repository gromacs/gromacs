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
    set(files
        utils.cpp
        load_table_f_2xnn.cpp
        load_table_f_v_2xnn.cpp
        load_table_f_4xn.cpp
        load_table_f_v_4xn.cpp
        arithmetic_one_input.cpp
        arithmetic_two_inputs.cpp
        arithmetic_three_inputs.cpp
        sincos.cpp
        logic.cpp
        anytrue.cpp
        blendzero.cpp
        masknotadd.cpp
        load.cpp
        set1.cpp
        setzero.cpp
        store.cpp
        )
    if (APPLE AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND
        C_COMPILER_VERSION AND C_COMPILER_VERSION VERSION_LESS 4.2)
            # A bug in Apple clang "4.1" prevents compiling two
            # templated forms of AlignedSimdAllocator in the same
            # compilation unit. I get the same error message as
            # https://code.google.com/p/address-sanitizer/issues/detail?id=132
            # but I do not understand why.

            # TODO I suspect the table loads are also broken because the
            # result SIMD vectors have to be allocated in aligned
            # memory??
        else()
        list(APPEND files
            load_interactions_2xnn.cpp
            load_interactions_4xn.cpp
            )
    endif()
    gmx_add_unit_test(SIMDTests simd-test
        ${files}
        )
endif()
