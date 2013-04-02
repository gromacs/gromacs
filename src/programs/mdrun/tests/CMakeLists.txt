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

#add_custom_target(mdrun-tests)

set(testname "mdrun-rerun-test")
set(testdir ${CMAKE_CURRENT_BINARY_DIR}/${testname})
file(MAKE_DIRECTORY ${testdir})

# input files
set(mdpfile ${CMAKE_CURRENT_SOURCE_DIR}/simple.mdp)
set(conffile ${CMAKE_CURRENT_SOURCE_DIR}/spc2.gro)
set(topfile ${CMAKE_CURRENT_SOURCE_DIR}/spc2.top)
set(rerunfile ${CMAKE_CURRENT_SOURCE_DIR}/spc2.trr)
# output files
set(tprfile ${testname}.tpr)

configure_file(mdrun-tests-main.cpp ${testname}/main.cpp)

set(MDRUN_SOURCES
    do_gct.c      gctio.c       genalg.c    ionize.c
    md.c          mdrun.c     membed.c
    pme_loadbal.c repl_ex.c     runner.c    xutils.c
    )
file(GLOB GROMPP_SOURCES ../../grompp/*.c)

unset(actual_sources)
foreach(mdrun_source_file ${MDRUN_SOURCES})
    list(APPEND actual_sources "../${mdrun_source_file}")
endforeach()
list(APPEND actual_sources ${GROMPP_SOURCES})

include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

add_executable(${testname} ${testname}/main.cpp ${actual_sources})
target_link_libraries(${testname} ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})

add_test(
    NAME ${testname}
    COMMAND ${testname}
    -s ${testname}
    -rerun ${rerunfile}
    -nt 1
    -ntomp 1
    WORKING_DIRECTORY ${testdir}
    )
set_tests_properties(${testname} PROPERTIES LABELS "Integration")

#add_dependencies(mdrun-tests ${testname}.tpr)
