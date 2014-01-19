#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
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

include_directories(${CMAKE_SOURCE_DIR}/src/gromacs/gmxpreprocess)

file(GLOB MDRUN_SOURCES mdrun/*.c mdrun/*.cpp)
# make an "object library" that we can re-use for multiple targets
add_library(mdrun_objlib OBJECT ${MDRUN_SOURCES})

if(GMX_FAHCORE)
    # The lack of a real source file here alongside the object library
    # may break some generators, according to CMake documentation. If
    # so, we can consider adding some dummy file to make it work.
    add_library(fahcore $<TARGET_OBJECTS:mdrun_objlib>)
elseif(GMX_BUILD_MDRUN_ONLY)
    add_executable(mdrun $<TARGET_OBJECTS:mdrun_objlib> mdrun_main.cpp)
    target_link_libraries(mdrun
        ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
    set_target_properties(mdrun PROPERTIES
        OUTPUT_NAME "mdrun${GMX_BINARY_SUFFIX}"
        COMPILE_FLAGS "${OpenMP_C_FLAGS}")
    install(TARGETS mdrun DESTINATION ${BIN_INSTALL_DIR} COMPONENT mdrun)
else()
    file(GLOB GMX_GENERAL_SOURCES gmx/*.c gmx/*.cpp)
    file(GLOB GMX_MAIN_SOURCES gmx/gmx.cpp gmx/legacymodules.cpp)
    list(REMOVE_ITEM GMX_GENERAL_SOURCES ${GMX_MAIN_SOURCES})
    if(GMX_X11)
        file(GLOB VIEW_SOURCES view/*.cpp)
    else()
        file(GLOB VIEW_SOURCES view/view.cpp)
    endif()
    add_library(gmx_objlib OBJECT ${GMX_GENERAL_SOURCES})
    add_library(view_objlib OBJECT ${VIEW_SOURCES})
    add_executable(gmx
        ${GMX_MAIN_SOURCES}
        $<TARGET_OBJECTS:gmx_objlib>
        $<TARGET_OBJECTS:mdrun_objlib>
        $<TARGET_OBJECTS:view_objlib>)
    target_link_libraries(gmx
        ${GMX_EXTRA_LIBRARIES} libgromacs ${GMX_EXE_LINKER_FLAGS})
    if(GMX_X11)
        target_link_libraries(gmx ${X11_LIBRARIES})
    endif()
    set_target_properties(gmx PROPERTIES
        OUTPUT_NAME "gmx${GMX_BINARY_SUFFIX}"
        COMPILE_FLAGS "${OpenMP_C_FLAGS}")

    install(TARGETS gmx
            RUNTIME DESTINATION ${BIN_INSTALL_DIR})
    if (GMX_SYMLINK_OLD_BINARY_NAMES)
        install(SCRIPT ${CREATE_LINKS_SCRIPT} COMPONENT links)
    endif ()

    if(BUILD_TESTING)
        add_subdirectory(mdrun/tests)
    endif(BUILD_TESTING)
endif()
