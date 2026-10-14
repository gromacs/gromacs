#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2021- The GROMACS Authors
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

if(GMX_HAVE_COLVARS)
    gmx_add_unit_test(ColvarsAppliedForcesUnitTest colvars_applied_forces-test
        CPP_SOURCE_FILES
            colvars.cpp
            colvarsoptions.cpp
            colvarspreprocessor.cpp
            colvarsforceprovider.cpp

            )
else()
    gmx_add_unit_test(ColvarsAppliedForcesUnitTest colvars_applied_forces-test
        CPP_SOURCE_FILES
            colvars.cpp
            colvarsoptions.cpp
        )
endif()

target_link_libraries(colvars_applied_forces-test PRIVATE
    applied_forces
    math
    topology
    mdlib
    mdtypes
    )
if (TARGET colvars_library)
    # Note that the design of Colvars requires that the GROMACS
    # MDModule include the definition of its types and perhaps define
    # some macros, so the GROMACS unit tests must depend on both the
    # module and transitively the headers and compile definitions used to
    # build it. Thus we obtain headers from the Colvars used to build
    # its MDModule
    target_include_directories(colvars_applied_forces-test PRIVATE
        $<TARGET_PROPERTY:colvars_library,INTERFACE_INCLUDE_DIRECTORIES>
        )
    target_compile_definitions(colvars_applied_forces-test PRIVATE
        $<TARGET_PROPERTY:colvars_library,COMPILE_DEFINITIONS>
        )
endif()
