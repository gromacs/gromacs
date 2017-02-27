#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2016, by the GROMACS development team, led by
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

# Helper functions for creating custom commands
#
# CMake semantics of add_custom_command() and add_custom_target() are not
# always very convenient or intuitive for creating commands that do not always
# run.  This file provides a gmx_add_custom_output_target() to simplify the
# task.  The function also provides some convenience features to remove code
# duplication in some parts of the GROMACS build system.
#
# Additionally, gmx_get_stamp_filename() and gmx_get_files_with_gitattribute()
# are provided independently to help in creating custom commands.

# Helper function to create a stamp file name for a target.
#
# Usage:
#   gmx_get_stamp_filename(<variable> <targetname>)
#
#   <variable>   - name of variable to receive the stamp name
#   <targetname> - name of target for which to generate the stamp name
#
# This is used internally by gmx_add_custom_target(... OUTPUT STAMP ...), but
# can also be called directly to create uniform stamp file names throughout the
# build system.
# <targetname> can be any string that is a part of a valid file name; it does
# not need to name an existing or to-be-created target.
function (gmx_get_stamp_filename variable targetname)
    set(_filename "${targetname}")
    if (NOT "${targetname}" MATCHES "stamp$")
        set(_filename "${targetname}-timestamp")
    endif()
    set(${variable} "${CMAKE_CURRENT_BINARY_DIR}/${_filename}.txt"
        PARENT_SCOPE)
endfunction()

# Helper function to tell gmx_add_custom_output_target() the name of an output
# file for a custom target.
#
# Usage:
#   gmx_set_custom_target_output(<targetname> <output>)
#
#   <targetname> - name of an existing custom target
#   <output>     - path to the output file produced by the custom target
#
# This is used internally by gmx_add_custom_output_target(), but can also be
# called for targets created with add_custom_target() to make them work with
# the dependency resolution mechanism in gmx_add_custom_output_target() if
# those targets for some reason cannot be created with that command.
function (gmx_set_custom_target_output targetname output)
    # Store the output file name in a custom property to be used in dependency
    # resolution later.
    if (NOT IS_ABSOLUTE ${output})
        set(output ${CMAKE_CURRENT_BINARY_DIR}/${output})
    endif()
    set_property(TARGET ${targetname}
        PROPERTY GMX_CUSTOM_TARGET_OUTPUT_FILE ${output})
endfunction()

# More flexible alternative to add_custom_command() and add_custom_target()
# for dependent custom commands.  It adds a few convenience features:
#   - Adds file-level dependencies between custom targets added with this
#     command such that if there is a target-level dependency, it also implies
#     that the custom command should always be run if the output file of the
#     dependency has been updated.
#   - Support for creating custom targets that produce stamp files whenever
#     they run successfully, so that other targets can depend on those stamp
#     files.
#
# Usage:
#   gmx_add_custom_output_target(<target> [RUN_ALWAYS] [ADD_FAST_TARGET]
#                                OUTPUT <STAMP | <output> >
#                                [COMMAND <command1> [<args1...>]]
#                                [COMMAND <command2> [<args2...>]]
#                                [WORKING_DIRECTORY <dir>]
#                                [DEPENDS <deps...>]
#                                [DEPENDS_FILE_LIST <list>]
#                                [COMMENT <comment>] [USES_TERMINAL])
#
#   <target>
#     - Name of the custom target to create.
#   RUN_ALWAYS
#     - Create the command such that it always runs, and may update the output
#       file in the process.
#       The dependencies listed with DEPENDS are ignored in this case.
#   ADD_FAST_TARGET
#     - In addition to creating <target>, create a secondary target
#       <target>-fast that always runs the same commands, but does not have
#       any dependencies.  Desired dependencies can be added separately using
#       add_dependencies().  This supports cases where some of the dependencies
#       are time-consuming to build, and it is possible to build the target
#       even without up-to-date dependencies for testing only that part of the
#       build.
#   OUTPUT
#     - Sets the name of the output file from this custom target.
#       Can be set to STAMP, in which case a stamp file name is automatically
#       generated and commands to touch that stamp whenever the target is made
#       are added.
#   COMMAND
#     - Passed to add_custom_command()/add_custom_target().
#       If OUTPUT STAMP is used, then a command to touch the stamp is
#       automatically added.  In this case, it is allowed to not specify any
#       commands explicitly.
#   WORKING_DIRECTORY
#     - Passed to add_custom_command()/add_custom_target()
#   COMMENT
#     - Passed to add_custom_command()/add_custom_target()
#   USES_TERMINAL
#     - Passed to add_custom_command()/add_custom_target()
#   DEPENDS
#     - Dependencies passed to add_custom_command().  Any targets in this list
#       that have been created with gmx_add_custom_output_target() are
#       automatically expanded such that the custom command always runs if the
#       output of the dependent target changes.  It is not necessary to list
#       those output files here explicitly.
#   DEPENDS_FILE_LIST
#     - Names of variables from which dependencies are added verbatim.
#       The target expansion described above is not performed for these
#       dependencies, and the value passed to the function is the name of a
#       list variable, not the list itself.  This provides much better
#       performance if the dependency list is a long list of source files.
#
# This function does not need a VERBATIM argument; that is always used when
# creating the underlying custom commands/targets.
function (gmx_add_custom_output_target targetname)
    # Parse the arguments
    # CMakeParseArguments is not suitable, since it does not support the use of
    # multiple COMMAND parameters that add_custom_target/command() supports.
    set(_option "")
    set(_command_args "")
    set(_deps "")
    set(_output "")
    set(_stamp OFF)
    set(_always OFF)
    set(_add_fast OFF)
    foreach (_arg ${ARGN})
        if ("x${_arg}" STREQUAL "xRUN_ALWAYS")
            set(_always ON)
        elseif ("x${_arg}" STREQUAL "xADD_FAST_TARGET")
            set(_add_fast ON)
        elseif ("x${_arg}" MATCHES "^x(OUTPUT|DEPENDS|DEPENDS_FILE_LIST)$")
            set(_option ${_arg})
        elseif ("x${_arg}" MATCHES "^x(COMMAND|COMMENT|WORKING_DIRECTORY|USES_TERMINAL)$")
            set(_option "PASS")
            list(APPEND _command_args "${_arg}")
        elseif ("x${_option}" STREQUAL "xDEPENDS")
            list(APPEND _deps "${_arg}")
            # If the dependency is a target created with this command, also add
            # the output file as a dependency.
            if (TARGET "${_arg}")
                get_property(_target_output
                    TARGET "${_arg}" PROPERTY GMX_CUSTOM_TARGET_OUTPUT_FILE)
                if (_target_output)
                    list(APPEND _deps ${_target_output})
                endif()
            endif()
        elseif ("x${_option}" STREQUAL "xPASS")
            list(APPEND _command_args "${_arg}")
        elseif ("x${_option}" STREQUAL "xDEPENDS_FILE_LIST")
            list(APPEND _deps ${${_arg}})
        elseif ("x${_option}" STREQUAL "xOUTPUT")
            if (_output)
                message(FATAL_ERROR "Multiple OUTPUTs not supported")
            endif()
            if ("x${_arg}" STREQUAL "xSTAMP")
                gmx_get_stamp_filename(_output ${targetname})
                set(_stamp ON)
            else()
                set(_output ${_arg})
            endif()
        else()
            message(FATAL_ERROR "Unknown option ${_arg}")
        endif()
    endforeach()
    # Add automatically a command to update the stamp if requested
    if (_stamp)
        list(APPEND _command_args COMMAND ${CMAKE_COMMAND} -E touch ${_output})
    endif()
    # Create the actual command as requested.
    if (NOT _always)
        add_custom_command(OUTPUT ${_output}
            ${_command_args} DEPENDS ${_deps} VERBATIM)
        add_custom_target(${targetname} DEPENDS ${_output})
    else()
        add_custom_target(${targetname} BYPRODUCTS ${_output}
            ${_command_args} VERBATIM)
    endif()
    # Store the output file name in a custom property to be used in dependency
    # resolution later.
    gmx_set_custom_target_output(${targetname} ${_output})
    # Create the fast target if requested.
    if (_add_fast)
        add_custom_target(${targetname}-fast ${_command_args} VERBATIM)
    endif()
endfunction()

# Gets list of files in the source tree with the given git attribute set
#
# Usage:
#   gmx_get_files_with_gitattribute(<variable> <attribute>)
#
#   <variable>  - name of variable to receive the file list
#   <attribute> - name of git attribute to use
#
# This is useful to generate a list of dependencies for custom commands when
# there is a long list of files that cannot be easily just globbed.
# Using git attributes allows keeping complicated logic out of the build
# system, as expressing such logic in a CMake script that would be reasonably
# fast to evaluate for a long file list is convenient or easy.
function (gmx_get_files_with_gitattribute variable attribute)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} ls-files
        COMMAND ${GIT_EXECUTABLE} check-attr ${attribute} --stdin
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE _files)
    string(REGEX MATCHALL "[^\n]*: ${attribute}: set\n"
           _files_with_attr "${_files}")
    string(REGEX REPLACE "([^;]*): ${attribute}: set\n" "${PROJECT_SOURCE_DIR}/\\1"
           _files "${_files_with_attr}")
    set(${variable} ${_files} PARENT_SCOPE)
endfunction()
