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
  file(STRINGS environment_script_output.txt output_lines)
  foreach(line IN LISTS output_lines)
    if (line MATCHES "^([a-zA-Z0-9_-]+)=(.*)$")
      set(ENV{${CMAKE_MATCH_1}} "${CMAKE_MATCH_2}")
    endif()
  endforeach()
endif()

file(TO_CMAKE_PATH "$ENV{GITHUB_WORKSPACE}" ccache_basedir)
set(ENV{CCACHE_BASEDIR} "${ccache_basedir}")
set(ENV{CCACHE_DIR} "${ccache_basedir}/.ccache")
set(ENV{CCACHE_COMPRESS} "true")
set(ENV{CCACHE_COMPRESSLEVEL} "6")
set(ENV{CCACHE_MAXSIZE} "600M")

execute_process(COMMAND ccache -p)
execute_process(COMMAND ccache -z)

execute_process(
  COMMAND cmake --build build
  RESULT_VARIABLE result-build
  OUTPUT_VARIABLE output-build
  ERROR_VARIABLE output-build
  ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
)
execute_process(
  COMMAND cmake --build build --target tests
  RESULT_VARIABLE result-build-test
  OUTPUT_VARIABLE output-build-test
  ERROR_VARIABLE output-build-test
  ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
)
if ((NOT result-build EQUAL 0) OR (NOT result-build-test EQUAL 0))
  string(REGEX MATCH "FAILED:.*$" error_message_build "${output-build}")
  string(REGEX MATCH "FAILED:.*$" error_message_build_test "${output-build-test}")
  string(REPLACE "\n" "%0A" error_message_build "${error_message_build}")
  string(REPLACE "\n" "%0A" error_message_build_test "${error_message_build_test}")
  message("::error::${error_message_build}")
  message("::error::${error_message_build_test}")
  message(FATAL_ERROR "Build failed")
endif()

execute_process(COMMAND ccache -s)

