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
install(FILES CMakeLists.txt.template
        DESTINATION ${DATA_INSTALL_DIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

file(GLOB_RECURSE GROMACS_HEADERS ${CMAKE_SOURCE_DIR}/include/*.h)
add_custom_command(OUTPUT gromacs
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/include gromacs 
    DEPENDS ${GROMACS_HEADERS})
add_custom_target(gromacs_include_links DEPENDS gromacs)

option(GMX_BUILD_TEMPLATE "Build gromacs template program" ON)
mark_as_advanced(GMX_BUILD_TEMPLATE)
# GMX_PREFER_STATIC_OPENMP=yes is a special case to build binaries
# to distribute and as the template is not installed it can be
# ignored.
# The template is build in a user-like environment, hence we use
# flags from PKG_CFLAGS. Again GMX_PREFER_STATIC_OPENMP=yes would
# need special link flags (OpenMP_LINKER_FLAGS), which are not
# very user-like.
if (GMX_BUILD_TEMPLATE AND NOT GMX_PREFER_STATIC_OPENMP)
    add_executable(template template.c)
    remove_definitions( -DHAVE_CONFIG_H )
    add_definitions("${PKG_CFLAGS}")
    target_link_libraries(template gmx)
    include_directories("${CMAKE_CURRENT_BINARY_DIR}")
    add_dependencies(template gromacs_include_links)
endif()

install(FILES README template.c Makefile.pkg
        DESTINATION ${DATA_INSTALL_DIR}/template
        COMPONENT development)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindGROMACS.cmakein
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake @ONLY)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/FindGROMACS.cmake
        DESTINATION ${DATA_INSTALL_DIR}/template/cmake
        COMPONENT development)

