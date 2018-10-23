#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2015,2016,2018, by the GROMACS development team, led by
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

include(CMakeParseArguments)

macro(gmx_init_sphinx_setup SPHINX_INPUT_DIR)
    set(_SPHINX_INPUT_ROOT ${SPHINX_INPUT_DIR})
    set(_SPHINX_INPUT_FILES)
    set(_SPHINX_IMAGE_CONVERSION_FILES)
endmacro()

macro(gmx_add_sphinx_input_file FILEPATH)
    list(APPEND _SPHINX_INPUT_FILES ${FILEPATH})
endmacro()

function(gmx_add_sphinx_source_files)
    set(_one_value_args FROM TO PREFIX)
    cmake_parse_arguments(ARG "" "${_one_value_args}" "FILES" ${ARGN})
    if (NOT ARG_FROM)
        set(ARG_FROM ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    if (ARG_TO)
        set(ARG_TO ${ARG_TO}/)
    endif()
    foreach(_file ${ARG_FILES})
        set(_source ${ARG_FROM}/${_file})
        get_filename_component(_filepath ${_file} DIRECTORY)
        get_filename_component(_filename ${_source} NAME)
        set(_targetdir ${_SPHINX_INPUT_ROOT}${ARG_TO}/${_filepath})
        set(_target ${_SPHINX_INPUT_ROOT}/${ARG_TO}${_file})
        if (NOT EXISTS ${_targetdir})
            file(MAKE_DIRECTORY ${_targetdir})
        endif()
        add_custom_command(
            OUTPUT  ${_target}
            COMMAND ${CMAKE_COMMAND} -E copy ${_source} ${_target}
            DEPENDS ${_source}
            COMMENT "Copying Sphinx input file ${ARG_PREFIX}${_file}"
            VERBATIM)
        list(APPEND _SPHINX_INPUT_FILES ${_target})
    endforeach()
    set(_SPHINX_INPUT_FILES "${_SPHINX_INPUT_FILES}" PARENT_SCOPE)
endfunction()

macro(gmx_remove_obsolete_sphinx_input_files IGNORE_PATTERN)
    file(GLOB_RECURSE _obsolete_sources ${_SPHINX_INPUT_ROOT}/*.rst)
    list(REMOVE_ITEM _obsolete_sources ${_SPHINX_INPUT_FILES})
    foreach(_file ${_obsolete_sources})
        file(RELATIVE_PATH _rel_path ${_SPHINX_INPUT_ROOT} ${_file})
        if (NOT _rel_path MATCHES "${IGNORE_PATTERN}")
            message(STATUS "Removing obsolete Sphinx input ${_rel_path}")
            file(REMOVE ${_file})
        endif()
    endforeach()
endmacro()

macro(gmx_add_sphinx_input_target TARGETNAME)
    gmx_add_custom_output_target(${TARGETNAME} OUTPUT STAMP
        DEPENDS ${_SPHINX_INPUT_FILES})
endmacro()

function(gmx_add_sphinx_image_conversion_files)
    set(_one_value_args FROM TO PREFIX)
    cmake_parse_arguments(ARG "" "${_one_value_args}" "FILES" ${ARGN})
    if (NOT ARG_FROM)
        set(ARG_FROM ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    if (ARG_TO)
        set(ARG_TO ${ARG_TO}/)
    endif()
    foreach(_file ${ARG_FILES})
        set(_source ${ARG_FROM}/${_file})
        get_filename_component(_filepath ${_file} DIRECTORY)
        get_filename_component(_filename ${_source} NAME)
        string(REGEX REPLACE "pdf" "png" _tmp ${_filename})
        set(_targetdir ${_SPHINX_INPUT_ROOT}/${ARG_TO}/${_filepath})
        set(_target ${_targetdir}/${_tmp})
        if (NOT EXISTS ${_targetdir})
            file(MAKE_DIRECTORY ${_targetdir})
        endif()
        add_custom_command(
            OUTPUT  ${_target}
            COMMAND ${ImageMagick_convert_EXECUTABLE} ${_source} -antialias -quality 03 -quiet -pointsize 12 -density 1200 -units PixelsPerInch  ${_target}
            DEPENDS ${_source}
            COMMENT "Converting Sphinx input graphics file ${_file} to png"
            VERBATIM)
        list(APPEND _SPHINX_IMAGE_CONVERSION_FILES ${_target})
    endforeach()
    set(_SPHINX_IMAGE_CONVERSION_FILES "${_SPHINX_IMAGE_CONVERSION_FILES}" PARENT_SCOPE)
endfunction()

macro(gmx_add_sphinx_image_conversion_target TARGETNAME)
    gmx_add_custom_output_target(${TARGETNAME} OUTPUT STAMP
        DEPENDS ${_SPHINX_IMAGE_CONVERSION_FILES})
endmacro()
