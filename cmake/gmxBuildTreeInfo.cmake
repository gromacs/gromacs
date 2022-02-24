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

# Retrieves information about the nature of the build tree
#
# The following variables are defined:
#   SOURCE_IS_SOURCE_DISTRIBUTION  The source tree is from a source tarball.
#   SOURCE_IS_GIT_REPOSITORY       The source tree is a git repository.
# Note that both can be false if the tree has been extracted, e.g., as a
# tarball directly from git.
# Additionally, the following variable is defined:
#   BUILD_IS_INSOURCE              The build is happening in-source.

#####################################################################
# Basic nature of the source tree

set(SOURCE_IS_GIT_REPOSITORY OFF)
set(SOURCE_IS_SOURCE_DISTRIBUTION OFF)
if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
    set(SOURCE_IS_GIT_REPOSITORY ON)
endif()
# This file is excluded from CPack source packages, but part of the repository,
# so it should get included everywhere else.
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/admin/.isreposource")
    set(SOURCE_IS_SOURCE_DISTRIBUTION ON)
endif()
set(BUILD_IS_INSOURCE OFF)
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    set(BUILD_IS_INSOURCE ON)
endif()

#####################################################################
# Location of other repositories (for development use only)

if (NOT DEFINED RELENG_PATH)
    if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/../releng")
        set(RELENG_PATH "${PROJECT_SOURCE_DIR}/../releng")
    endif()
endif()
set(RELENG_PATH "${RELENG_PATH}" CACHE PATH "Path to releng repository")
mark_as_advanced(RELENG_PATH)
