#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2020- The GROMACS Authors
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

# The legacy installed API consists of headers that are not considered
# to be maintainable in terms of a stable API specification. These headers
# will no longer be available to install in a future release.

# Note: Any usage requirements that should be transitive should be added to
# this INTERFACE target. Compiler and linker options (that do not need to
# be propagated when linking to the `common` target) can be added directly
# to the `common` target.
add_library(legacy_api INTERFACE)
target_include_directories(legacy_api INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                           $<INSTALL_INTERFACE:include>)
if(BUILD_SHARED_LIBS AND WIN32)
  target_compile_definitions(legacy_api INTERFACE
                             $<BUILD_INTERFACE:libgromacs_EXPORTS>)
endif()

configure_file(version.h.cmakein include/gromacs/version.h)

# Allow an export target for legacy_api since it is a dependency of the
# installable libgromacs target.
install(TARGETS legacy_api
        EXPORT libgromacs
        INCLUDES DESTINATION include)

if(GMX_INSTALL_LEGACY_API)
    # Install public header directories.
    install(DIRECTORY include/gromacs
            DESTINATION include)

    # Install "configured" files from the build tree.
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/gromacs/version.h
            DESTINATION include/gromacs)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/gromacs/libgromacs_export.h
            DESTINATION include/gromacs)
endif()
