#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#

include_directories(${OpenMM_INCLUDE_DIR})
link_directories(${OpenMM_LIBRARY_DIR}) 
# with this define no evn.var. is needed with OPENMM_PLUGIN_DIR
# if the same OpenMM installation is used for running and building 
add_definitions( -DOPENMM_PLUGIN_DIR="${OpenMM_PLUGIN_DIR}" ) 
file(TO_CMAKE_PATH ${OpenMM_PLUGIN_DIR} _path)
add_library(openmm_api_wrapper STATIC ${CMAKE_SOURCE_DIR}/src/contrib/openmm_wrapper.cpp)
target_link_libraries(openmm_api_wrapper ${OpenMM_LIBRARIES})

list(REMOVE_ITEM MDRUN_SOURCES mdrun.c runner.c)
list(APPEND MDRUN_SOURCES
    ${CMAKE_SOURCE_DIR}/src/contrib/md_openmm.c
    ${CMAKE_SOURCE_DIR}/src/contrib/mdrun_openmm.c
    ${CMAKE_SOURCE_DIR}/src/contrib/runner_openmm.c
    )

# this is to circumvent the following MSVC error: 
# warning LNK4098: defaultlib 'LIBCMT' conflicts with use of other libs
# fatal error LNK1169: one or more multiply defined symbols found
if(GMX_OPENMM AND MSVC)
    set_target_properties(mdrun PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src/gmxlib/gpu_utils)

set_source_files_properties(main.c PROPERTIES LANGUAGE CXX)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set_source_files_properties(main.c PROPERTIES COMPILE_FLAGS "-x c++")
endif()
