#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2019, by the GROMACS development team, led by
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

file(GLOB MODULARSIMULATOR_SOURCES *.cpp)

add_library(modularsimulator OBJECT ${MODULARSIMULATOR_SOURCES})
gmx_target_compile_options(modularsimulator)
target_compile_definitions(modularsimulator PRIVATE HAVE_CONFIG_H)
target_include_directories(modularsimulator SYSTEM BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/src/external/thread_mpi/include)

if(GMX_OPENMP)
    # Explicitly set properties for modular simulator module to compile with openmp
    set_target_properties(modularsimulator PROPERTIES COMPILE_OPTIONS $<TARGET_PROPERTY:OpenMP::OpenMP_CXX,INTERFACE_COMPILE_OPTIONS>)
endif()

if (WIN32)
    gmx_target_warning_suppression(modularsimulator /wd4244 HAS_NO_MSVC_LOSSY_CONVERSION)
    gmx_target_warning_suppression(modularsimulator /wd4996 HAS_NO_MSVC_UNSAFE_FUNCTION)
else()
    # Several std::move uses are redundant in C++14, but clang before
    # 3.9 had a bug that needed them. gcc 9 can warn about their use,
    # which we need to supress. This suppression should be removed
    # when GROMACS requires clang 3.9 or higher.
    gmx_target_warning_suppression(modularsimulator "-Wno-redundant-move" HAS_NO_REDUNDANT_MOVE)
endif()

list(APPEND libgromacs_object_library_dependencies modularsimulator)
set(libgromacs_object_library_dependencies ${libgromacs_object_library_dependencies} PARENT_SCOPE)
