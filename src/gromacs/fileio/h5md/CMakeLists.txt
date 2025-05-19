#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2025- The GROMACS Authors
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

# Set up the module library
add_library(h5md OBJECT)

gmx_manage_hdf5()

if(GMX_USE_HDF5)
    target_sources(h5md PRIVATE
        h5md.cpp
        h5md_attribute.cpp
        h5md_dataset.cpp
        h5md_error.cpp
        h5md_group.cpp
        h5md_guard.cpp
        h5md_read.cpp
        h5md_write.cpp
    )

    target_include_directories(h5md SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(h5md PUBLIC ${HDF5_LIBRARIES})
else()
    # We do not need to compile H5md implementation files if we aren't building with it
    target_sources(h5md PRIVATE h5md.cpp)
endif()

target_link_libraries(h5md PRIVATE
    common
    legacy_api
    legacy_modules
    math
    utility)

if (BUILD_SHARED_LIBS)
    set_target_properties(h5md PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Set required compile flags
gmx_target_compile_options(h5md)
target_compile_definitions(h5md PRIVATE HAVE_CONFIG_H)

if (BUILD_TESTING)
     add_subdirectory(tests)
endif()
