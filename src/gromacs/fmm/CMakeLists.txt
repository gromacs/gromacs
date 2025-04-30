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

add_library(fmm OBJECT)

gmx_manage_fmm()

target_sources(fmm PRIVATE
    fmm_mdmodule.cpp
    fmm_mdpoptions.cpp
    fmmoptions.cpp
)

if(GMX_USE_EXT_FMM)
    target_sources(fmm PRIVATE
        fmmforceprovider.cpp
    )

else()
    target_sources(fmm PRIVATE
        fmmforceprovider_stub.cpp
    )

endif()

if (BUILD_SHARED_LIBS)
    set_target_properties(fmm PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

gmx_target_compile_options(fmm)
target_compile_definitions(fmm PRIVATE HAVE_CONFIG_H)

target_link_libraries(fmm PRIVATE
    common
    legacy_modules
    math
    mdtypes
    utility
    )

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
