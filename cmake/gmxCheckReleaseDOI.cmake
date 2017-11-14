#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2017, by the GROMACS development team, led by
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

# Include newly registered DOI into Gromacs for new releases.

# Configure script to include the DOI registered during Jenkins
# release build into GROMACS source code
# Checks if the values for GMX_SOURCE_DOI and GMX_MANUAL_DOI
# are set during build time. If so, then the file gmxDOIVersion.cmake.cmakein
# gets configures with those values to be included into the source
# tarball. If not, and if the configured file exists (only in tarballs
# originating from release builds), the strings get included into 
# the subsequent builds. If the file does not exists, empty strings 
# will be used that generate message that it is not a release build.

#
# Paul Bauer (paul.bauer.q@gmail.com)

# master variables controlling build and submission
option(GMX_GET_RELEASE_DOI "Activate for a build of the release tarball containing doi from Zenodo" OFF)
mark_as_advanced(GMX_GET_RELEASE_DOI)

# set file path variables so they are available everywhere below
set(INFILE "${CMAKE_SOURCE_DIR}/cmake/gmxDOIVersion.cmake.cmakein")
set(OUTFILE "${CMAKE_SOURCE_DIR}/cmake/gmxDOIVersion.cmake")

# Only do this if we are building the release binary
# Instead, check for existing gmxDOIVersion.cmake file and use it to populate
# the variables for source and manual, if build from the source code tarball
# if the file does not exist, the build is either from a git repository or
# an archive not generated from a release build. In this case, the file
# is generated with empty variables so that the version stays marked
# as a build not originating from a release
if((NOT GMX_GET_RELEASE_DOI) AND (NOT (EXISTS ${OUTFILE})))
    # check if the populated gmxDOIVersion file exists or not
    # if not, we are not originating from a release build and the 
    # doi strings should be empty
    set(GMX_MANUAL_DOI "" CACHE INTERNAL "reserved doi for GROMACS manual")
    set(GMX_SOURCE_DOI "${GMX_MANUAL_DOI}" CACHE INTERNAL "reserved doi for GROMACS manual")
    configure_file(${INFILE} ${OUTFILE})
else()
    if(GMX_GET_RELEASE_DOI)
        # this build was started from the workflow that also registers the doi
        # so we need to have the values for both source and manual doi passed
        # to CMake to configure the gmxDOIVersion.cmake.cmakein file with them later
        if(NOT (DEFINED GMX_SOURCE_DOI) OR NOT (DEFINED GMX_MANUAL_DOI))
            # both strings need to be set at least with dummy variables for this to work
            MESSAGE(FATAL_ERROR "GMX_SOURCE_DOI or GMX_MANUAL_DOI not set during release build")
        endif()
        configure_file(${INFILE} ${OUTFILE})
    endif()
    # if the previous statement was not true, we started from a release
    # build archive with the gmxDOIVersion.cmake file generated there
    # or have already build a non release version previously, so we keep the empty variables
endif()

# include the new file after it became populated
include(gmxDOIVersion)
