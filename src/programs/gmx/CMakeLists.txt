#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2012, by the GROMACS development team, led by
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

add_executable(gmx gmx.cpp legacymodules.cpp gmx_main.cpp)
target_link_libraries(gmx libgromacs gmxana ${GMX_EXE_LINKER_FLAGS})
set_target_properties(gmx PROPERTIES OUTPUT_NAME "gmx${GMX_BINARY_SUFFIX}")

configure_file(CreateLinks.cmake.cmakein CreateLinks.cmake @ONLY)
set(CREATE_LINKS_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/CreateLinks.cmake)
# TODO: Currently this rule is always run and recreates the symlinks;
# it would be best to just adapt the rest of the build such that this is
# never needed (for now, this makes the binary directory seem like it was
# earlier).
add_custom_target(links ALL
    ${CMAKE_COMMAND}
        -DBINARY_DIRECTORY=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        -DQUIETLY=TRUE
        -P ${CREATE_LINKS_SCRIPT}
    DEPENDS gmx
    COMMENT "Creating symlinks for alternative binary names"
    VERBATIM)

# TODO: Create man pages (currently not working since the old mechanism
# requires a target for each binary to which one can add a post-build role, and
# the symlinks do not have such targets).
# May be worth to wait until #685 and #969 are resolved; after these, it should
# be possible to handle everything uniformly instead of having several
# different CMake-side mechanisms.

install(TARGETS gmx
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})
if (GMX_SYMLINK_OLD_BINARY_NAMES)
    install(SCRIPT ${CREATE_LINKS_SCRIPT} COMPONENT links)
endif ()

if (BUILD_TESTING)
    add_subdirectory(legacytests)
endif (BUILD_TESTING)
