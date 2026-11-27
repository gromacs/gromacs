#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
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


# Add Lepton library, which is developed and distributed as part of OpenMM.
# See src/external/lepton/README.md for details


gmx_option_multichoice(GMX_USE_LEPTON
    "Build the Lepton library interfaced with GROMACS"
    INTERNAL
    INTERNAL NONE)
mark_as_advanced(GMX_USE_LEPTON)

function(gmx_manage_lepton)
  if(GMX_USE_LEPTON STREQUAL "INTERNAL")

    set(LEPTON_DIR "${CMAKE_SOURCE_DIR}/src/external/lepton")

    file(GLOB LEPTON_SOURCES ${LEPTON_DIR}/src/*.cpp)
    add_library(lepton_objlib OBJECT ${LEPTON_SOURCES})
    target_include_directories(lepton_objlib SYSTEM PRIVATE ${LEPTON_DIR}/include)

    # For simplicity, Lepton is built as a static library for now
    target_compile_options(lepton_objlib PRIVATE -DLEPTON_BUILDING_STATIC_LIBRARY)
    set_target_properties(lepton_objlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

    add_library(lepton INTERFACE)
    target_link_libraries(lepton INTERFACE
      $<BUILD_INTERFACE:lepton_objlib>
      $<TARGET_OBJECTS:lepton_objlib>
      )
    target_include_directories(lepton SYSTEM INTERFACE $<BUILD_INTERFACE:${LEPTON_DIR}/include>)

  else()

    # Add a dummy target
    add_library(lepton INTERFACE)
  endif()
endfunction()

if(GMX_USE_LEPTON STREQUAL "INTERNAL")
    set(GMX_HAVE_LEPTON 1 CACHE INTERNAL "Is Lepton found?")
else()
    set(GMX_HAVE_LEPTON 0 CACHE INTERNAL "Is Lepton found?")
endif()
mark_as_advanced(GMX_HAVE_LEPTON)
