#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2012,2013,2014,2015,2018,2019, by the GROMACS development team, led by
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

file(GLOB SELECTION_SOURCES *.cpp)
file(GLOB SCANNER_SOURCES scanner.cpp parser.cpp)
list(REMOVE_ITEM SELECTION_SOURCES ${SCANNER_SOURCES})
# Add the non-generated sources to libgromacs, so that they have only the normal warning suppressions
target_sources(libgromacs PRIVATE ${SELECTION_SOURCES})

add_library(scanner OBJECT ${SCANNER_SOURCES})
gmx_target_compile_options(scanner)
target_compile_definitions(scanner PRIVATE HAVE_CONFIG_H)
# Compilers give lots of warnings on these generated files.  These
# are awkward to suppress inline, so we do it in the compilation
# command (after testing that the compiler supports the
# suppressions).
if (WIN32)
    gmx_target_warning_suppression(scanner /wd4101 HAS_NO_MSVC_UNUSED)
    gmx_target_warning_suppression(scanner /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION)
    gmx_target_warning_suppression(scanner /wd4267 HAS_NO_MSVC_LOSSY_CONVERSION_SIZE_T_TO_INT)
else()
    gmx_target_warning_suppression(scanner -Wno-unused HAS_NO_UNUSED)
    gmx_target_warning_suppression(scanner -Wno-unused-parameter HAS_NO_UNUSED_PARAMETER)
    gmx_target_warning_suppression(scanner -Wno-missing-declarations HAS_NO_MISSING_DECLARATIONS)
    gmx_target_warning_suppression(scanner -wd1419 HAS_DECL_IN_SOURCE)
endif()
list(APPEND libgromacs_object_library_dependencies scanner)
set(libgromacs_object_library_dependencies ${libgromacs_object_library_dependencies} PARENT_SCOPE)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
