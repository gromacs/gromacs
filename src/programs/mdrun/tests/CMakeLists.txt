#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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

set(testname "MdrunTests")
set(exename "mdrun-tests")

set(mdrun_dir "${CMAKE_CURRENT_SOURCE_DIR}/..")
file(GLOB MDRUN_SOURCES ${mdrun_dir}/*.c ${mdrun_dir}/*.cpp)
file(GLOB MDRUN_SOURCES ../*.c ../*.cpp)
set(gmx_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../gmx")
file(GLOB GMX_EXECUTABLE_SOURCES ${gmx_dir}/grompp.c ${gmx_dir}/convparm.c ${gmx_dir}/tomorse.c)
list(APPEND actual_sources ${MDRUN_SOURCES})
list(APPEND actual_sources ${GMX_EXECUTABLE_SOURCES})

include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

# TODO By design, multiple tools might be run when using MPI, so
# finalization with gmx::finalize() should be avoided until the end of
# all tests in an executable. The only way to ensure it is called is
# to use a different main() function from the unit tests. We could do
# that, but can probably cope without finalization, too.
gmx_build_unit_test(
    ${testname}
    ${exename}
    # files with code for tests
    rerun.cpp
    # files with code for test fixtures
    moduletest.cpp
    # files for code for grompp and mdrun
    ${actual_sources}
    )
gmx_register_integration_test(
    ${testname}
    ${exename}
    )
