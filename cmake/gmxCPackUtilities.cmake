#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2015, by the GROMACS development team, led by
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

# Helper functions to encapsulate and modularize usage of CPack
#
# This file is intended to be the only place that directly sets CPack variables
# (in gmx_cpack_write_config()), and other parts of the build system should
# only call the functions declared here to set up the packaging.

# Initialize the machinery that collects CPack information during other
# build system generation
#
# This function should be called before other functions from this file.
function (gmx_cpack_init)
    # Add the source tree to the source package.
    # If we would not set CPACK_SOURCE_INSTALLED_DIRECTORIES, this is what
    # CPack would set as default.
    set_property(GLOBAL PROPERTY GMX_CPACK_SOURCE_INSTALLED_DIRECTORIES
        ${PROJECT_SOURCE_DIR} /)
endfunction()

# Add a generated directory to be included in the source package.
#
# Usage:
#   gmx_cpack_add_generated_source_directory(<dir> [DESTINATION <dest>])
#
#   <dir>   Name of directory to include.
#           Relative paths are interpreted relative to current build dir.
#   <dest>  Path in the source package where files from <dir> will be put.
#           If not set, the files are put in the corresponding location in
#           the source tree.
#
# By default, CPack source archives includes all files from the source tree.
# This function adds a directory from the build tree to be packaged into the
# source archive.  These are used for content GROMACS generates as part of the
# configuration or build.
# The values end up in CPACK_SOURCE_INSTALLED_DIRECTORIES, which is a list of
# pairs of names of source and destination directories.
function (gmx_cpack_add_generated_source_directory DIR)
    include(CMakeParseArguments)
    set(_one_value_args DESTINATION)
    cmake_parse_arguments(ARG "" "${_one_value_args}" "" ${ARGN})
    if (ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "Unknown arguments: ${ARG_UNPARSED_ARGUMENTS}")
    endif()
    set(_dir ${DIR})
    if (NOT IS_ABSOLUTE ${_dir})
        set(_dir ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
    endif()
    if (ARG_DESTINATION)
        set(_dest ${ARG_DESTINATION})
    else()
        file(RELATIVE_PATH _dest ${PROJECT_BINARY_DIR} ${_dir})
    endif()
    set_property(GLOBAL APPEND PROPERTY GMX_CPACK_SOURCE_INSTALLED_DIRECTORIES
        ${_dir} ${_dest})
endfunction()

# Write CPack configuration files
#
# This function should be called at the end of the main CMakeLists.txt, after
# all other calls to functions in this file.
# CPack also automatically populates the list of components based on components
# used in installation rules, so it should come after all install() commands.
function (gmx_cpack_write_config)
    # Set basic package information.
    set(CPACK_PACKAGE_NAME    "gromacs")
    set(CPACK_PACKAGE_VENDOR  "gromacs.org")
    set(CPACK_PACKAGE_CONTACT "gmx-users@gromacs.org")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
        "GROMACS - a toolkit for high-performance molecular simulation")
    # Set version info.
    set(CPACK_PACKAGE_VERSION_MAJOR ${GMX_VERSION_MAJOR})
    set(CPACK_PACKAGE_VERSION_MINOR ${GMX_VERSION_MINOR})
    set(CPACK_PACKAGE_VERSION_PATCH ${GMX_VERSION_PATCH})
    set(CPACK_PACKAGE_VERSION       ${GMX_VERSION_STRING})
    # Add various text resources for some installers.
    set(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/admin/InstallWelcome.txt")
    # Its GPL/LGPL, so they do not have to agree to a license for mere usage,
    # but some installers require this...
    set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
    set(CPACK_RESOURCE_FILE_README  "${PROJECT_SOURCE_DIR}/admin/InstallInfo.txt")

    # Our custom config file that is run by CPack for each generator, used to
    # check for prerequisites of the packaging.
    set(CPACK_PROJECT_CONFIG_FILE "${PROJECT_SOURCE_DIR}/CPackInit.cmake")

    # Settings specific to source packages.
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_IGNORE_FILES
        "\\\\.isreposource$;\\\\.git/;\\\\.gitignore$;\\\\.gitattributes;")
    # Get the list of directories added with gmx_cpack_add_generated_source_directory()
    get_property(CPACK_SOURCE_INSTALLED_DIRECTORIES
        GLOBAL PROPERTY GMX_CPACK_SOURCE_INSTALLED_DIRECTORIES)

    # Propagate additional values for CPackInit.cmake to use.
    # CPack includes all variables starting with CPACK_ into the generated
    # config files that are included by CPack.
    set(CPACK_GMX_BUILD_HELP "${GMX_BUILD_HELP}")

    # Generate the CPack configuration files.
    include(CPack)
endfunction()
