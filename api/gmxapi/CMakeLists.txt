#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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

# This list file provides the Gromacs::gmxapi cmake module.

##########################
# Set up public interface.

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
set(GMXAPI_MINOR 3)
set(GMXAPI_PATCH 0)
set(GMXAPI_RELEASE ${GMXAPI_MAJOR}.${GMXAPI_MINOR}.${GMXAPI_PATCH})

add_library(gmxapi)
set_target_properties(gmxapi PROPERTIES
                      VERSION_MAJOR ${GMXAPI_MAJOR}
                      VERSION_MINOR ${GMXAPI_MINOR}
                      VERSION_PATCH ${GMXAPI_PATCH}
                      RELEASE ${GMXAPI_RELEASE})


if (GMX_LIB_MPI)
    # GROMACS is built against an MPI library.
    set(_gmx_mpi_type "library")
    target_link_libraries(gmxapi PRIVATE MPI::MPI_CXX)
elseif(GMX_THREAD_MPI)
    # GROMACS is built with its internal thread-MPI implementation.
    set(_gmx_mpi_type "tmpi")
else()
    # GROMACS is built without MPI or thread-MPI.
    set(_gmx_mpi_type "none")
endif ()
define_property(TARGET PROPERTY MPI
                BRIEF_DOCS "MPI capability of the GROMACS library."
                FULL_DOCS "Values of 'library', 'tmpi', or 'none' indicate the configure-time options of the GROMACS
                 library build.")
set_target_properties(gmxapi PROPERTIES
                      MPI ${_gmx_mpi_type})
unset(_gmx_mpi_type)

##########################
# Define public interface.
#
# The include directory should be mostly empty so that we can use it internally as
# the public interface include directory during build and testing.
configure_file(include/gmxapiversion.h.in include/gmxapi/version.h)
configure_file(include/resourceassignment.cmakein.h include/gmxapi/mpi/resourceassignment.h)
target_include_directories(gmxapi PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                           $<INSTALL_INTERFACE:include>
                           )

add_subdirectory(cpp)

###############################
# Install the public interface.
#
# If any item begins in a generator expression it must evaluate to a full path,
# so we can't just use something like $<TARGET_PROPERTIES:gmxapiPublicHeaders,SOURCES>.
# Instead, we use a canonical list defined in the parent scope.

# Install public header directories.
install(DIRECTORY include/gmxapi
        DESTINATION include)

# Install "configured" files from the build tree.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/gmxapi/version.h
        DESTINATION include/gmxapi)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/gmxapi/mpi/resourceassignment.h
        DESTINATION include/gmxapi/mpi)

