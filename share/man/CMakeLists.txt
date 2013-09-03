#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
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

MACRO(TODAY RESULT)
    IF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%F" OUTPUT_VARIABLE ${RESULT}
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    ELSE()
        set(${RESULT} "????-??-??")
    ENDIF()
ENDMACRO(TODAY)

if(NOT SOURCE_IS_SOURCE_DISTRIBUTION)
    if(NOT CMAKE_CROSSCOMPILING)
        # Defaults to off until there is a reasonable solution for cases
        # where the binary does not execute on the compilation host.
        option(GMX_BUILD_MANPAGES "Build man pages" OFF)
    else()
        message(STATUS "Building the man pages is not available when "
            "cross-compiling the developer version from git")
    endif()
else()
    #make sure source package contains all man pages
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/man1/gmx-view.1")
        message(FATAL_ERROR "Man pages are missing from source package.")
    endif()
endif()
mark_as_advanced(GMX_BUILD_MANPAGES)

if(GMX_BUILD_MANPAGES)
    # create gromacs.7
    FILE(READ "${CMAKE_SOURCE_DIR}/admin/programs.txt" contents)

    # Convert file contents into a CMake list. First escape ;
    STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
    STRING(REGEX REPLACE "\n" ";" contents "${contents}")

    set(PROGMANPAGES "")
    foreach(line ${contents})
        if (${line} MATCHES "^HEAD\\|")
            string(REGEX REPLACE "^HEAD\\|" "" DESC ${line})
            set(PROGMANPAGES "${PROGMANPAGES}.Sh \"${DESC}\"\n.IX Subsection \"${DESC}\"\n.Vb\n.ta 16n\n")
        elseif(${line} MATCHES "^END$")
            set(PROGMANPAGES "${PROGMANPAGES}.Ve\n")
        elseif(${line} MATCHES "\\|")
            string(REGEX REPLACE "\\|" "\t" line ${line})
            set(PROGMANPAGES "${PROGMANPAGES}\\&  ${line}\n")
        else()
            message(WARNING "Incorrectly formated line \"${line}\" in programs.txt")
        endif()
    endforeach()
    TODAY(TODAYS_DATE)
    configure_file(man7/gromacs.7.cmakein man7/gromacs.7)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/man7/gromacs.7
        DESTINATION ${MAN_INSTALL_DIR}/man7
        COMPONENT man)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/man1)
    # TODO: It could be nicer to not do this as part of the ALL target,
    # but instead on demand as part of make install.
    add_custom_target(man ALL
        gmx -quiet help -export man
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating man pages for programs"
        VERBATIM)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/man1
        DESTINATION ${MAN_INSTALL_DIR}
        COMPONENT man)
    set_directory_properties(PROPERTIES
        ADDITIONAL_MAKE_CLEAN_FILES man1)

# man pages are only available if they are either build or this is a source archive
elseif(SOURCE_IS_SOURCE_DISTRIBUTION)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/man7/gromacs.7
        DESTINATION ${MAN_INSTALL_DIR}/man7
        COMPONENT man)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/man1
        DESTINATION ${MAN_INSTALL_DIR}
        COMPONENT man)
endif()
