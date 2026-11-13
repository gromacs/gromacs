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

if (GMX_USE_HDF5)
    set(H5MD_TEST_SOURCE_FILES
        h5md.cpp
        h5mdtestbase.cpp
        h5md_attribute.cpp
        h5md_compound_type.cpp
        h5md_datasetbase.cpp
        h5md_datasetbuilder.cpp
        h5md_fixeddataset.cpp
        h5md_framedataset.cpp
        h5md_framedatasetbuilder.cpp
        h5md_group.cpp
        h5md_guard.cpp
        h5md_indexgroup.cpp
        h5md_particleblock.cpp
        h5md_sz3_compression.cpp
        h5md_timedatablock.cpp
        h5md_topologyutils.cpp
        h5md_type.cpp
        h5md_util.cpp
        binary_information.cpp
    )
else()
    set(H5MD_TEST_SOURCE_FILES h5md_disabled.cpp)
endif()

gmx_add_unit_test(H5mdTests h5md-test CPP_SOURCE_FILES ${H5MD_TEST_SOURCE_FILES})
target_link_libraries(h5md-test PRIVATE
    h5md
    topology
)
if (GMX_HAVE_SZ3)
    target_link_libraries(h5md-test PRIVATE hdf5sz3)
endif()

# Set whether we are building with SZ3 as a compile-time definition in target source files
set_source_files_properties(binary_information.cpp h5md_sz3_compression.cpp
    PROPERTIES
    COMPILE_DEFINITIONS "GMX_HAVE_SZ3=$<BOOL:${GMX_HAVE_SZ3}>"
)
