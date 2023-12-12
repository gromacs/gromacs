#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2014- The GROMACS Authors
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

gmx_add_gtest_executable(gmxpreprocess-test
    CPP_SOURCE_FILES
        editconf.cpp
        genconf.cpp
        genion.cpp
        genrestr.cpp
        gen_maxwell_velocities.cpp
        gpp_atomtype.cpp
        gpp_bond_atomtype.cpp
        grompp_directives.cpp
        insert_molecules.cpp
        massrepartitioning.cpp
        readir.cpp
        solvate.cpp
        topdirs.cpp
        )
gmx_register_gtest_test(GmxPreprocessTests gmxpreprocess-test SLOW_TEST)
target_link_libraries(gmxpreprocess-test PRIVATE
        gmxpreprocess
        topology
        )

# Currently these can be slow to run in CI, so they are in
# several test binaries.

set(exename pdb2gmx1-test)
gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
       pdb2gmx.cpp
       )
target_compile_definitions(${exename} PRIVATE OPLSAA=1 GROMOS=0 AMBER=0 CHARMM=0)
gmx_register_gtest_test(Pdb2gmx1Test ${exename} SLOW_TEST)

set(exename pdb2gmx2-test)
gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        pdb2gmx.cpp
        )
target_compile_definitions(${exename} PRIVATE OPLSAA=0 GROMOS=1 AMBER=0 CHARMM=0)
gmx_register_gtest_test(Pdb2gmx2Test ${exename} SLOW_TEST)

set(exename pdb2gmx3-test)
gmx_add_gtest_executable(${exename}
    CPP_SOURCE_FILES
        pdb2gmx.cpp
        )
target_compile_definitions(${exename} PRIVATE OPLSAA=0 GROMOS=0 AMBER=1 CHARMM=1)
gmx_register_gtest_test(Pdb2gmx3Test ${exename} SLOW_TEST)
