#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018, by the GROMACS development team, led by
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

# Build test input files to test directory.
# Note that CMake interprets `gmx` in `COMMAND gmx` as the executable produced
# by the `add_executable` command that defined the `gmx` target.

set(GMXAPI_TEST_TIMESTEP 0.001953125) # 2e-9
# Note that the CMake configure_file below substitutes GMXAPI_TEST_TIMESTEP for the mdp `dt` value.
# \todo replace command line wrapper with API utilities or test framework facilities.
# Directly manage input data instead of relying on the following in-repo input files.
# \todo When migrating to an API-based testing configuration, use GMXAPI_TEST_TIMESTEP for `dt`
configure_file(data/spc-and-methanol.mdp.in
               ${CMAKE_CURRENT_BINARY_DIR}/data/spc-and-methanol.mdp)
add_custom_command(
    OUTPUT topol.tpr
    COMMAND gmx grompp
        # Note: at time of writing, these are insufficient to suppress all output.
        -quiet
        -nocopyright
        -nov
        -o topol.tpr
        -f ${CMAKE_CURRENT_BINARY_DIR}/data/spc-and-methanol.mdp
        -c ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.gro
        -n ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.ndx
        -p ${CMAKE_CURRENT_SOURCE_DIR}/data/spc-and-methanol.top
        # The following syntax is platform dependent to redirect stdout and stderr
        # to the null device, but initial gmxapi support is limited, anyway.
        # TODO: Find a quieter way to generate this input data on all platforms.
        > /dev/null 2>&1
    DEPENDS gmx
            data/spc-and-methanol.gro
            ${CMAKE_CURRENT_BINARY_DIR}/data/spc-and-methanol.mdp
            data/spc-and-methanol.top
            data/spc-and-methanol.ndx
)

add_custom_target(spc_and_methanol_tpr
                  DEPENDS topol.tpr
                  COMMENT "Preparing gmxapi test reference data."
                  )

configure_file(testingconfiguration.in.cpp testingconfiguration.cpp)

add_library(gmxapi-test-support SHARED testingconfiguration.cpp)
add_dependencies(gmxapi-test-support spc_and_methanol_tpr)
target_include_directories(gmxapi-test-support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
