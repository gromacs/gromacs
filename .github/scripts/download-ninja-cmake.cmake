#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2022- The GROMACS Authors
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

set(cmake_version $ENV{CMAKE_VERSION})
set(ninja_version $ENV{NINJA_VERSION})

message(STATUS "Using host CMake version: ${CMAKE_VERSION}")
message(STATUS "Using RUNNER_OS: $ENV{RUNNER_OS}")
message(STATUS "Using RUNNER_ARCH: $ENV{RUNNER_ARCH}")

if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
    set(ninja_suffix "win.zip")
    set(cmake_suffix "win64-x64.zip")
    set(cmake_dir "cmake-${cmake_version}-win64-x64/bin")
elseif ("$ENV{RUNNER_OS}" STREQUAL "macOS")
    set(ninja_suffix "mac.zip")
    set(cmake_suffix "Darwin-x86_64.tar.gz")
    set(cmake_dir "cmake-${cmake_version}-Darwin-x86_64/CMake.app/Contents/bin")
elseif ("$ENV{RUNNER_OS}" STREQUAL "Linux" AND "$ENV{RUNNER_ARCH}" STREQUAL "ARM64")
    set(ninja_suffix "linux.zip")
    set(cmake_suffix "linux-aarch64.tar.gz")
    set(cmake_dir "cmake-${cmake_version}-linux-aarch64/bin")
endif()

set(ninja_url "https://github.com/ninja-build/ninja/releases/download/v${ninja_version}/ninja-${ninja_suffix}")
file(DOWNLOAD "${ninja_url}" ./ninja.zip SHOW_PROGRESS)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ./ninja.zip)

set(cmake_url "https://github.com/Kitware/CMake/releases/download/v${cmake_version}/cmake-${cmake_version}-${cmake_suffix}")
file(DOWNLOAD "${cmake_url}" ./cmake.zip SHOW_PROGRESS)
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xvf ./cmake.zip)

# Add to PATH environment variable
file(TO_CMAKE_PATH "$ENV{GITHUB_WORKSPACE}/${cmake_dir}" cmake_dir)
set(path_separator ":")
if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
    set(path_separator ";")
endif()
file(APPEND "$ENV{GITHUB_PATH}" "$ENV{GITHUB_WORKSPACE}${path_separator}${cmake_dir}")

if (NOT "$ENV{RUNNER_OS}" STREQUAL "Windows")
    execute_process(
    COMMAND chmod +x ninja
    COMMAND chmod +x ${cmake_dir}/cmake
    )
endif()
