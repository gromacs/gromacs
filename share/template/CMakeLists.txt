#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2011- The GROMACS Authors
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

add_executable(template template.cpp)
if (WIN32)
    gmx_target_warning_suppression(template /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION_DOUBLE_TO_REAL)
    gmx_target_warning_suppression(template /wd4267 HAS_NO_MSVC_LOSSY_CONVERSION_SIZE_T_TO_INT)
endif()
# This should be removable once object libraries can directly use target_link_libraries
# with CMake 3.12, #3290
target_include_directories(template SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/src/external)
target_link_libraries(template PRIVATE libgromacs legacy_modules legacy_api ${GMX_EXE_LINKER_FLAGS})

set(DOCUMENTATION_HTTP_URL_BASE
    https://manual.gromacs.org/nightly/index.html)
if (SOURCE_IS_SOURCE_DISTRIBUTION)
    set(DOCUMENTATION_HTTP_URL_BASE
        http://manual.gromacs.org/documentation/${GMX_VERSION_STRING})
endif()
configure_file(README.cmakein README @ONLY)

install(FILES CMakeLists.txt.template
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template
        RENAME CMakeLists.txt
        COMPONENT development)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/README template.cpp Makefile.pkg
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template
        COMPONENT development)

install(FILES cmake/FindGROMACS.cmake
        DESTINATION ${GMX_INSTALL_GMXDATADIR}/template/cmake
        COMPONENT development)
