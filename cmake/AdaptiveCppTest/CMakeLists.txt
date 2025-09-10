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

cmake_minimum_required(VERSION 3.28)
cmake_policy(VERSION 3.28)

project(AdaptiveCppTest)
find_package(adaptivecpp QUIET)
if(NOT adaptivecpp_FOUND)
    find_package(hipsycl REQUIRED)
endif()

# This project will fail compilation if it can't use AdaptiveCpp compiler,
# and will emit a warning for each of the known compilation targets:
# - "GMX_ACPP_TEST_HAVE_CUDA_TARGET" for CUDA
# - "GMX_ACPP_TEST_HAVE_HIP_TARGET" for HIP
# - "GMX_ACPP_TEST_HAVE_LEVELZERO_TARGET" for deprecated LevelZero
# - "GMX_ACPP_TEST_HAVE_GENERIC_TARGET" for generic/SSCP

add_executable(AdaptiveCppTest main.cpp)
add_sycl_to_target(TARGET AdaptiveCppTest SOURCES main.cpp)
