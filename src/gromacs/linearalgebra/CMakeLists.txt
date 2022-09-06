#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
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

add_library(linearalgebra OBJECT)
target_sources(linearalgebra PRIVATE
        eigensolver.cpp
        gmx_arpack.cpp
        matrix.cpp
        nrjac.cpp
        sparsematrix.cpp
)

if (NOT GMX_EXTERNAL_BLAS)
    file(GLOB BLAS_SOURCES gmx_blas/*.cpp)
    target_sources(linearalgebra PRIVATE ${BLAS_SOURCES})
endif()

if (NOT GMX_EXTERNAL_LAPACK)
    file(GLOB LAPACK_SOURCES gmx_lapack/*.cpp)
    target_sources(linearalgebra PRIVATE ${LAPACK_SOURCES})
endif()

if (BUILD_SHARED_LIBS)
    set_target_properties(linearalgebra PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(linearalgebra PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories(linearalgebra PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

gmx_target_compile_options(linearalgebra)
target_compile_definitions(linearalgebra PRIVATE HAVE_CONFIG_H)
# The linearalgebra code is all considered external, and we will
# not keep it free of warnings. Any compiler suppressions required
# should be added here.
if (WIN32)
    gmx_target_warning_suppression(linearalgebra /wd4305 HAS_NO_MSVC_INITIALIZING_FLOAT_FROM_DOUBLE)
    gmx_target_warning_suppression(linearalgebra /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION_DOUBLE_TO_REAL)
    gmx_target_warning_suppression(linearalgebra /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
else()
    # We usually enable this warning, but here it is being used
    # correctly to construct strings for use with an API that does
    # not expect null termination of C strings.
    gmx_target_warning_suppression(linearalgebra -Wno-stringop-truncation HAS_NO_STRINGOP_TRUNCATION)
endif()
target_link_libraries(linearalgebra PRIVATE legacy_api)
target_link_libraries(linearalgebra PRIVATE common)
# TODO: Remove `legacy_modules` once specific modules are explicitly linked.
target_link_libraries(linearalgebra PRIVATE legacy_modules utility)
