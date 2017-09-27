#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

gmx_add_unit_test(SimdUnitTests simd-test
                  bootstrap_loadstore.cpp
		  base.cpp
                  data.cpp
                  simd.cpp
		  simd_floatingpoint.cpp
                  simd_floatingpoint_util.cpp
                  simd_vector_operations.cpp
                  simd_math.cpp
		  simd_integer.cpp
                  simd4.cpp
                  simd4_floatingpoint.cpp
                  simd4_vector_operations.cpp
                  simd4_math.cpp
                  scalar.cpp
                  scalar_util.cpp
                  scalar_math.cpp)

# Add tests for expressions which are supposed to not compile when the
# build configuration supports a real implementation of the SIMD
# module. See tests/simd_ambiguous.cpp for documentation.
set(AMBIGUOUS_COMPILATION_EXPECTED_TO_FAIL TRUE)
if(GMX_SIMD_ACTIVE STREQUAL "NONE")
    set(AMBIGUOUS_COMPILATION_EXPECTED_TO_FAIL FALSE)
endif()
foreach(TEST_PREC float double)
    foreach(TEST_FUNC exp exp2 log inv cos sin sqrt)
        set(TEST_NAME simd_ambiguous_${TEST_PREC}_${TEST_FUNC})
        add_executable(${TEST_NAME} simd_ambiguous.cpp)
        set_target_properties(${TEST_NAME} PROPERTIES
                              EXCLUDE_FROM_ALL TRUE
                              EXCLUDE_FROM_DEFAULT_BUILD TRUE)
        target_compile_definitions(${TEST_NAME} PRIVATE TEST_FUNC=${TEST_FUNC} TEST_PREC=${TEST_PREC})
        add_test(NAME ${TEST_NAME}
            COMMAND ${CMAKE_COMMAND} --build . --target ${TEST_NAME} --config $<CONFIGURATION>
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL ${AMBIGUOUS_COMPILATION_EXPECTED_TO_FAIL})
    endforeach()
endforeach()
