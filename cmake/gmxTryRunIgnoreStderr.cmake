#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

# - A (highly) simplified version of try_run() that allows us to compile
#   a single source file with user-provided additional definitions,
#   but in contrast to the standard CMake version this one only captures
#   stdout and ignores stderr (since some braindead compilers appear to think
#   it's a good idea to echo unsolicited output to stderr when running a binary).
#
# gmx_try_run_ignore_stderr(RUN_RESULT_VAR
#                           COMPILE_RESULT_VAR
#                           <bindir>
#                           <srcfile>
#                           [COMPILE_DEFINITIONS <flags>]
#                           [COMPILE_OUTPUT_VARIABLE comp]
#                           [RUN_OUTPUT_VARIABLE run]
#                           [ARGS <arg1> <arg2> ...]
#
# The available options have similar roles as the CMake-provided try_run(), but
# be aware that this version might be a lot more fragile.
#

include(CMakeParseArguments)

function(gmx_try_run_ignore_stderr RUN_RESULT_VAR COMPILE_RESULT_VAR BINDIR SRCFILE)

    if(NOT DEFINED ${RUN_RESULT_VAR})

        set(_options ) # empty
        set(_one_value_args COMPILE_DEFINITIONS COMPILE_OUTPUT_VARIABLE RUN_OUTPUT_VARIABLE)
        set(_multi_value_args ARGS)
        cmake_parse_arguments(ARG "${_options}" "${_one_value_args}" "${_multi_value_args}" ${ARGN})

        if(NOT DEFINED ${COMPILE_RESULT_VAR})
            # First compile, capture stdout but ignore stderr.
            # Save output binary with a name that includes COMPILE_RESULT_VAR, so we can cache it.
            try_compile(${COMPILE_RESULT_VAR}
                        "${BINDIR}"
                        "${SRCFILE}"
                        COMPILE_DEFINITIONS "${ARG_COMPILE_DEFINITIONS}"
                        OUTPUT_VARIABLE COMPILE_OUTPUT
                        COPY_FILE ${CMAKE_BINARY_DIR}/CMakeFiles/TryRunNoStderr_${COMPILE_RESULT_VAR}${CMAKE_EXECUTABLE_SUFFIX})

            if(ARG_COMPILE_OUTPUT_VARIABLE)
                set(${ARG_COMPILE_OUTPUT_VARIABLE} "${COMPILE_OUTPUT}" PARENT_SCOPE)
            endif()
        endif(NOT DEFINED ${COMPILE_RESULT_VAR})

        # If it compiled, try to execute it, but ignore stderr.
        # This might use the cached output binary.
        if(${${COMPILE_RESULT_VAR}})
            execute_process(COMMAND ${CMAKE_BINARY_DIR}/CMakeFiles/TryRunNoStderr_${COMPILE_RESULT_VAR}${CMAKE_EXECUTABLE_SUFFIX} ${ARG_ARGS}
                            RESULT_VARIABLE ${RUN_RESULT_VAR}
                            OUTPUT_VARIABLE RUN_OUTPUT
                            ERROR_QUIET)

            if(ARG_RUN_OUTPUT_VARIABLE)
                set(${ARG_RUN_OUTPUT_VARIABLE} "${RUN_OUTPUT}" PARENT_SCOPE)
            endif()

            set(${RUN_RESULT_VAR} ${${RUN_RESULT_VAR}} PARENT_SCOPE)

        endif(${${COMPILE_RESULT_VAR}})

    endif(NOT DEFINED ${RUN_RESULT_VAR})

endfunction()
