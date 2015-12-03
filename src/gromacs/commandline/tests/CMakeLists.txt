#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

set(EXECUTABLE_EXTENSION "")
if (GMX_NATIVE_WINDOWS OR GMX_CYGWIN)
    set(EXECUTABLE_EXTENSION ".exe")
endif()
set(PATH_SEARCH_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/test-bin)
file(MAKE_DIRECTORY ${PATH_SEARCH_TEST_DIR}/bin)
file(WRITE ${PATH_SEARCH_TEST_DIR}/bin/test-exe${EXECUTABLE_EXTENSION}
     "Test executable for path searching")
if (UNIX)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            test-exe
            ${PATH_SEARCH_TEST_DIR}/bin/test-rel-link)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${PATH_SEARCH_TEST_DIR}/bin/test-exe
            ${PATH_SEARCH_TEST_DIR}/bin/test-abs-link)
endif()

gmx_add_unit_test(CommandLineUnitTests commandline-test
                  cmdlinehelpmodule.cpp
                  cmdlinehelpwriter.cpp
                  cmdlinemodulemanager.cpp
                  cmdlinemodulemanagertest.cpp
                  cmdlineparser.cpp
                  cmdlineprogramcontext.cpp
                  pargs.cpp
                  $<TARGET_OBJECTS:onlinehelp-test-shared>)

if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
    # This suppression stops a very verbose cascade of messages about the
    # mocks, which is probably a compiler issue.
    #   1540-2924 (W) Cannot pass an argument of non-POD class type "const gmx::CommandLineHelpContext" through ellipsis.
    set_property(SOURCE cmdlinehelpmodule.cpp cmdlinemodulemanager.cpp cmdlinemodulemanagertest.cpp
                 PROPERTY COMPILE_FLAGS "-qsuppress=1540-2924")
endif()
