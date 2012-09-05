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
                -DEXENAME=${EXENAME}
                -P ${CMAKE_SOURCE_DIR}/cmake/CreateManPage.cmake)
        install(FILES ${MAN1_PATH}/${EXENAME}.1 DESTINATION 
            ${MAN_INSTALL_DIR}/man1 OPTIONAL)
    elseif(NOT EXISTS "${CMAKE_SOURCE_DIR}/admin/.isreposource")
        install(FILES ${CMAKE_SOURCE_DIR}/man/man1/${EXENAME}.1 DESTINATION 
            ${MAN_INSTALL_DIR}/man1)
    endif()
endfunction ()
