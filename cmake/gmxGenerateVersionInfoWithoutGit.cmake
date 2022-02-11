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

# Generate GROMACS release build version information, as well as build
# directory integrity checking.

# This script generates version information for a build from a release tarball
# source tree based on the release information in the released tarball.
# It is assumed that by default the script is run in cmake script mode.
# If *not* called in script mode but used in generating cache variables,
# GEN_VERSION_INFO_INTERNAL has to be set ON.
#
# The following variables have to be previously defined:
# PROJECT_VERSION               - hard-coded version string (generated info is appended)
# PROJECT_SOURCE_DIR            - top level source directory
# VERSION_CMAKEIN               - path to an input template file
# VERSION_OUT                   - path to the output file
#
# Output:
# VERSION_OUT is configured from the input VERSION_CMAKEIN
# using the variables listed below.
#
# GMX_VERSION_STRING_FULL       - version string
#
# Paul Bauer (paul.bauer.q@gmail.com)
# (authors of git Version of the script that this is based on below)
# Szilard Pall (pszilard@cbr.su.se)
# Teemu Murtola (teemu.murtola@gmail.com)

# Check input variables.
if("${PROJECT_VERSION}" STREQUAL "")
    message(FATAL_ERROR "PROJECT_VERSION undefined!")
endif()
if (NOT EXISTS "${PROJECT_SOURCE_DIR}")
    message(FATAL_ERROR "Project source directory ${PROJECT_SOURCE_DIR} does not exist")
endif()
if ("${VERSION_CMAKEIN}" STREQUAL "")
    message(FATAL_ERROR "Missing input parameter VERSION_CMAKEIN!")
endif()
if ("${VERSION_OUT}" STREQUAL "")
    message(FATAL_ERROR "Missing input parameter VERSION_OUT!")
endif()

# Prepare version string to populate
set(GMX_VERSION_STRING_FULL          ${PROJECT_VERSION})

# Change path depending on if this is a source distribution (e.g. release tarball)
# or just a source directory that is managed by something else, like an IDE.
# If the git executable isn't found by CMake, there will not be version info even
# if the .git folder is present and SOURCE_IS_GIT_REPOSITORY is true.
# Don't issue warnings in this case.
set(GMX_VERSION_STRING_FULL "${GMX_VERSION_STRING_FULL}")

# Generate the output file.
configure_file(${VERSION_CMAKEIN} ${VERSION_OUT})
