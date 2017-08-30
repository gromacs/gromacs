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

set(GMXAPI_USERDOC_OUTPUT ${CMAKE_BINARY_DIR}/docs/html/doxygen/api-user)
set(GMXAPI_DEVDOC_OUTPUT ${CMAKE_BINARY_DIR}/docs/html/doxygen/api-dev)
set(GMXAPI_WEBDOC_OUTPUT ${CMAKE_BINARY_DIR}/docs/html/gmxapi)

file(MAKE_DIRECTORY
     ${GMXAPI_USERDOC_OUTPUT}
     ${GMXAPI_DEVDOC_OUTPUT}
     ${GMXAPI_WEBDOC_OUTPUT})

get_target_property(GMXAPI_DEV gmxapi SOURCE_DIR)
string (REPLACE ";" " " GMXAPI_PUBLIC "${GMXAPI_PUBLIC_HEADERS}")

configure_file(Doxyfile.in Doxyfile @ONLY)
configure_file(Doxyfile-dev.in Doxyfile-dev @ONLY)

# todo: the docs builds should probably be moved to separate targets under the appropriate subdirectories.

################################
# C++ external API documentation
################################
find_package(Doxygen)
if (DOXYGEN_FOUND)
	add_custom_target( gmxapi_cppdocs
	    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
	    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	    COMMENT "Generating API user documentation for ${GMXAPI_USERDOC_OUTPUT}"
	    VERBATIM )
	add_custom_target( gmxapi_cppdocs_dev
	    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile-dev
	    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	    COMMENT "Generating API developer documentation for ${GMXAPI_DEVDOC_OUTPUT}"
	    VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
