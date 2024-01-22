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

if ("$ENV{RUNNER_OS}" STREQUAL "Windows" AND NOT "x$ENV{ENVIRONMENT_SCRIPT}" STREQUAL "x")
  execute_process(
    COMMAND "$ENV{ENVIRONMENT_SCRIPT}" && set
    OUTPUT_FILE environment_script_output.txt
  )
  file(STRINGS environment_script_output.txt output_lines)
  foreach(line IN LISTS output_lines)
    if (line MATCHES "^([a-zA-Z0-9_-]+)=(.*)$")
      set(ENV{${CMAKE_MATCH_1}} "${CMAKE_MATCH_2}")
    endif()
  endforeach()
endif()

set(path_separator ":")
if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
  set(path_separator ";")
endif()
if ("$ENV{RUNNER_OS}" STREQUAL "macOS")
  # macOS supports RDTSCP, but there are issues with GitHub runners, #4896
  set(RDTSCP_VAR OFF)
else()
  set(RDTSCP_VAR ON)
endif()
set(ENV{PATH} "$ENV{GITHUB_WORKSPACE}${path_separator}$ENV{PATH}")

message(STATUS "Using GPU_VAR: $ENV{GPU_VAR}")

execute_process(
  COMMAND cmake
    -S .
    -B build
    -D CMAKE_BUILD_TYPE=$ENV{BUILD_TYPE}
    -G Ninja
    -D CMAKE_MAKE_PROGRAM=ninja
    -D CMAKE_C_COMPILER_LAUNCHER=ccache
    -D CMAKE_CXX_COMPILER_LAUNCHER=ccache
    -D GMX_COMPILER_WARNINGS=ON
    -D GMX_DEFAULT_SUFFIX=OFF
    -D GMX_GPU=$ENV{GPU_VAR}
    -D GMX_SIMD=None
    -D GMX_FFT_LIBRARY=FFTPACK
    -D GMX_OPENMP=$ENV{OPENMP_VAR}
    -D REGRESSIONTEST_DOWNLOAD=ON
    -D GMX_USE_RDTSCP=${RDTSCP_VAR}
  RESULT_VARIABLE result
)
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Bad exit status")
endif()
