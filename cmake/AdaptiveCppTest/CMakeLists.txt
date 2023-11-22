#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2021- The GROMACS Authors
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

cmake_minimum_required(VERSION 3.18)
project(AdaptiveCppTest)
find_package(adaptivecpp QUIET)
if(NOT adaptivecpp_FOUND)
    find_package(hipsycl REQUIRED)
endif()

# By default, this project will fail compilation if it can't use AdaptiveCpp/hipSYCL compiler.
# By defining CHECK_x_TARGET we can also make it fail if the specified backend is not among the compilation targets.

if (CHECK_CUDA_TARGET)
    add_compile_definitions(CHECK_CUDA_TARGET)
endif()
if (CHECK_HIP_TARGET)
    add_compile_definitions(CHECK_HIP_TARGET)
endif()
if (CHECK_LEVELZERO_TARGET)
    add_compile_definitions(CHECK_LEVELZERO_TARGET)
endif()
if (CHECK_GENERIC_TARGET)
    add_compile_definitions(CHECK_GENERIC_TARGET)
endif()
add_executable(AdaptiveCppTest main.cpp)
add_sycl_to_target(TARGET AdaptiveCppTest SOURCES main.cpp)
