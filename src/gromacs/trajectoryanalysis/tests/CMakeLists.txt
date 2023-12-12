#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2010- The GROMACS Authors
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

gmx_add_gtest_executable(trajectoryanalysis-test
    CPP_SOURCE_FILES
        # Infrastructure
        moduletest.cpp
        # Tests
        angle.cpp
        clustsize.cpp
        cmdlinerunner.cpp
        convert_trj.cpp
        distance.cpp
        dssp.cpp
        extract_cluster.cpp
        freevolume.cpp
        gyrate.cpp
        hbond.cpp
        msd.cpp
        pairdist.cpp
        rdf.cpp
        sasa.cpp
        select.cpp
        scattering.cpp
        surfacearea.cpp
        topologyinformation.cpp
        trajectory.cpp
        unionfind.cpp
        $<TARGET_OBJECTS:mdrun_objlib>
        )
gmx_register_gtest_test(TrajectoryAnalysisUnitTests trajectoryanalysis-test SLOW_TEST)
target_link_libraries(trajectoryanalysis-test PRIVATE analysisdata-test-shared)
target_link_libraries(trajectoryanalysis-test PRIVATE common)

add_executable(test_selection ${UNITTEST_TARGET_OPTIONS} test_selection.cpp)
target_link_libraries(test_selection PRIVATE
                      common
                      libgromacs
                      ${GMX_EXE_LINKER_FLAGS})
# TODO: Explicitly link to module dependencies.
target_link_libraries(test_selection PRIVATE legacy_modules utility)
