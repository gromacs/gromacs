#  This file is part of Gromacs        Copyright (c) 1991-2008
#  David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.

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
