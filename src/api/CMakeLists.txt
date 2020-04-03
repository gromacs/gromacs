#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

# Note: GROMACS releases have a single-integer monotonic version in GMX_API_VERSION
# and LIBRARY_VERSION annotates the shared object for libgromacs. GMXAPI versioning
# is not synchronized to releases and may increment faster or slower.
#
# Prior to 0.1, GMXAPI patch levels are used to mark short term development cycles
# and allow compatibility checks for client software of the early releases.
#
# gmxapi 0.2 will be the first release candidate for gmxapi 1.0 and will attempt
# to establish compatibility guarantees consistent with semantic versioning.
# (https://semver.org). When the API is deemed suitably stable, gmxapi 1.0 should
# be tagged. Official GROMACS releases should be mappable to a distinct gmxapi
# release string. For roadmap details, see https://gitlab.com/gromacs/gromacs/-/issues/2585
set(GMXAPI_MAJOR 0)
set(GMXAPI_MINOR 1)
set(GMXAPI_PATCH 0)
set(GMXAPI_RELEASE ${GMXAPI_MAJOR}.${GMXAPI_MINOR}.${GMXAPI_PATCH})

###############################
# New public C++ API components
###############################

# Define a list of the public headers to be installed and documented. Use
# absolute paths in the build tree. Path component before `gmxapi` can be
# stripped as appropriate by consumers of this list.
set(GMXAPI_PUBLIC_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/context.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/exceptions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/gmxapi.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/gmxapicompat.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/gromacsfwd.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/md.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/session.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/status.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/system.h
    ${CMAKE_CURRENT_BINARY_DIR}/cpp/include/gmxapi/version.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/compat/mdparams.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/compat/tpr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/gmxapi/md/mdmodule.h
    )

add_subdirectory(cpp)

add_subdirectory(docs)
