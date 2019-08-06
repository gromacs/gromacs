#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2018,2019, by the GROMACS development team, led by
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

file(GLOB LINEARALGEBRA_SOURCES *.cpp)

if (NOT GMX_EXTERNAL_BLAS)
    file(GLOB BLAS_SOURCES gmx_blas/*.cpp)
endif()

if (NOT GMX_EXTERNAL_LAPACK)
    file(GLOB LAPACK_SOURCES gmx_lapack/*.cpp)
endif()

set(LINEARALGEBRA_SOURCES
    ${LINEARALGEBRA_SOURCES} ${BLAS_SOURCES} ${LAPACK_SOURCES})

add_library(linearalgebra OBJECT ${LINEARALGEBRA_SOURCES})
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
list(APPEND libgromacs_object_library_dependencies linearalgebra)
set(libgromacs_object_library_dependencies ${libgromacs_object_library_dependencies} PARENT_SCOPE)
