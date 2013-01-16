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
#
MACRO(TODAY RESULT)
    IF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%F" OUTPUT_VARIABLE ${RESULT} OUTPUT_STRIP_TRAILING_WHITESPACE)
    ELSE()
        set(${RESULT} "????-??-??")
    ENDIF()
ENDMACRO(TODAY)

if(GMX_BUILD_MANPAGES)
    set(MAN1_PATH ${CMAKE_BINARY_DIR}/man/man1)
    file(MAKE_DIRECTORY ${MAN1_PATH})

    #create gromacs.7
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
    configure_file(${CMAKE_SOURCE_DIR}/man/man7/gromacs.7.cmakein ${CMAKE_BINARY_DIR}/man/man7/gromacs.7)
    install(FILES ${CMAKE_BINARY_DIR}/man/man7/gromacs.7 DESTINATION
        ${MAN_INSTALL_DIR}/man7)
#man-pages are only avalaible if they are either build or this is a source archive
elseif(NOT EXISTS "${CMAKE_SOURCE_DIR}/admin/.isreposource")
    install(FILES ${CMAKE_SOURCE_DIR}/man/man7/gromacs.7 DESTINATION
        ${MAN_INSTALL_DIR}/man7)
endif()

function (gmx_add_man_page EXENAME)
    if(GMX_BUILD_MANPAGES)
        file(STRINGS ${CMAKE_SOURCE_DIR}/admin/programs.txt DESC 
            REGEX "^${EXENAME}\\|" LIMIT_COUNT 1)
        #Regex breaks with a "|" in description. Cmake doesn't support 
        #non-greedy regex.
        string(REGEX REPLACE "^.*\\|" "" DESC "${DESC}")
        if(DESC STREQUAL "")
            message(WARNING "Missing description for ${EXENAME}")
        endif()
        add_custom_command(TARGET ${EXENAME} POST_BUILD 
            #The redirect is a hack to avoid showing copyright. 
            #Ideally -quiet would also cause programs to not print copyright.
            COMMAND ${CMAKE_COMMAND} -DINFILE=${EXENAME}${GMX_BINARY_SUFFIX}.nroff 
                -DOUTFILE=${MAN1_PATH}/${EXENAME}.1 -DDESC=" - ${DESC}"
                -DEXENAME=${EXENAME}${GMX_BINARY_SUFFIX}
                -P ${CMAKE_SOURCE_DIR}/cmake/CreateManPage.cmake)
        install(FILES ${MAN1_PATH}/${EXENAME}.1 DESTINATION 
            ${MAN_INSTALL_DIR}/man1 OPTIONAL)
    elseif(NOT EXISTS "${CMAKE_SOURCE_DIR}/admin/.isreposource")
        install(FILES ${CMAKE_SOURCE_DIR}/man/man1/${EXENAME}.1 DESTINATION 
            ${MAN_INSTALL_DIR}/man1)
    endif()
endfunction ()
