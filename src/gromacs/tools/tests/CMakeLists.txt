#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2018- The GROMACS Authors
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

gmx_add_gtest_executable(tool-test
    CPP_SOURCE_FILES
        dump.cpp
        helpwriting.cpp
        make_ndx.cpp
        report_methods.cpp
        trjconv.cpp
        convert-tpr.cpp
        )
gmx_register_gtest_test(ToolUnitTests tool-test SLOW_TEST)
target_link_libraries(tool-test PRIVATE math)

gmx_add_gtest_executable(tool-test-with-leaks
        CPP_SOURCE_FILES
        convert-tpr-with-leaks.cpp
        )
gmx_register_gtest_test(ToolWithLeaksUnitTests tool-test-with-leaks SLOW_TEST IGNORE_LEAKS)
target_link_libraries(tool-test-with-leaks PRIVATE math)
