#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013, by the GROMACS development team, led by
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

# The module defines the following variables:
#   VMD_EXECUTABLE - path to vmd command
#   GMX_VMD_PLUGIN_PATH - path to vmd plugins

message(STATUS "Checking for suitable VMD version")
find_program(VMD_EXECUTABLE NAMES vmd PATH_SUFFIXES bin
    DOC "VMD command")

#set search path in increasing priority:
# default path, vmd binary path, enviroment variable
set(VMD_PATHS "/usr/local/lib/vmd/plugins/*/molfile/")
if(VMD_EXECUTABLE)
    file(STRINGS "${VMD_EXECUTABLE}" VMDDIR REGEX "^defaultvmddir=.*$")
    string(REGEX REPLACE "(^.*=\"?|\"$)" "" VMDDIR "${VMDDIR}")
    list(INSERT VMD_PATHS 0 "${VMDDIR}/plugins/*/molfile/")
endif()
if(NOT "$ENV{VMDDIR}" STREQUAL "")
    list(INSERT VMD_PATHS 0 "$ENV{VMDDIR}/plugins/*/molfile/")
endif()

#xyz is just an example. Any other molfile plugin could be used.
#But some require extra link flags. VMD uses ".so" even on Windows.
find_library(VMDXYZPLUGIN NAME "xyzplugin.so"
    PATHS ${VMD_PATHS})

if (VMDXYZPLUGIN)
    try_run(TESTVMD TESTVMD_COMPILED ${CMAKE_BINARY_DIR}
        "${CMAKE_SOURCE_DIR}/cmake/TestVMD.c"
        CMAKE_FLAGS "-DLINK_LIBRARIES=${CMAKE_DL_LIBS}"
            "-DINCLUDE_DIRECTORIES=${CMAKE_SOURCE_DIR}/src/gmxlib"
        COMPILE_DEFINITIONS "-DGMX_USE_PLUGINS"
        RUN_OUTPUT_VARIABLE TESTVMD_OUTPUT ARGS ${VMDXYZPLUGIN})
endif()

if(NOT TESTVMD EQUAL 0)
    if (NOT VMDXYZPLUGIN)
        message(STATUS "VMD plugins not found. Path to VMD can be set with VMDDIR.")
    elseif(NOT TESTVMD_COMPILED)
        message(STATUS "Could not compile VMD version check")
    elseif(TESTVMD EQUAL 1)
        message(STATUS "Could not load VMD plugin ${VMDXYZPLUGIN}: ${TESTVMD_OUTPUT}")
    elseif(TESTVMD EQUAL 5)
        message(STATUS "VMD plugin ${VMDXYZPLUGIN} too old. VMD 1.8.6 required.")
    else()
        message(STATUS "Could not identify VMD version of ${VMDXYZPLUGIN}. Error: ${TESTVMD}")
    endif()
    # This permits GROMACS to avoid hard-coding a fall-back path.
    # Fall-back is useful in case VMD is installed later.
    set(GMX_VMD_PLUGIN_PATH "/usr/local/lib/vmd/plugins/*/molfile"
        CACHE PATH "Path to VMD plugins for molfile I/O" FORCE)
else()
    get_filename_component(VMD_PLUGIN_PATH ${VMDXYZPLUGIN} PATH)
    message(STATUS "VMD version of ${VMD_PLUGIN_PATH} is suitable")
    set(GMX_VMD_PLUGIN_PATH ${VMD_PLUGIN_PATH}
        CACHE PATH "Path to VMD plugins for molfile I/O" FORCE)
endif()
mark_as_advanced(GMX_VMD_PLUGIN_PATH)
#Nothing is rerun unless GMX_VMD_PLUGIN_PATH is set to NO. Clean-up all.
unset(VMDXYZPLUGIN CACHE)
unset(VMD_EXECUTABLE CACHE)
unset(TESTVMD CACHE)
unset(TESTVMD_COMPILED CACHE)
unset(VMD_PATHS)
unset(VMD_PLUGIN_PATH)
unset(VMDDIR)
