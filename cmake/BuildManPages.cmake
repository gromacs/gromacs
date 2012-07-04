if(GMX_BUILD_MANPAGES)
    set(MAN1_PATH ${CMAKE_BINARY_DIR}/man/man1)
    file(MAKE_DIRECTORY ${MAN1_PATH})
endif()
function (gmx_add_man_page EXENAME)
    if(GMX_BUILD_MANPAGES)
        file(STRINGS ${CMAKE_SOURCE_DIR}/admin/programs.txt DESC 
            REGEX "^${EXENAME}\\|")
        #Regex breaks with a "|" in description. Cmake doesn't support 
        #non-greedy regex.
        string(REGEX REPLACE "^.*\\|" "" DESC "${DESC}")
        if(DESC STREQUAL "")
            message(WARNING "Missing description for ${EXENAME}")
        endif()
        add_custom_command(TARGET ${EXENAME} POST_BUILD 
            #The redirect is a hack to avoid showing copyright. 
            #Ideally -quiet would also cause programs to not print copyright.
            COMMAND ${EXENAME} -quiet -man nroff 2>${EXENAME}.err
            COMMAND ${CMAKE_COMMAND} -DINFILE=${EXENAME}${GMX_BINARY_SUFFIX}.nroff 
                -DOUTFILE=${MAN1_PATH}/${EXENAME}.1 -DDESC=" - ${DESC}"
                -P ${CMAKE_SOURCE_DIR}/cmake/Filter.cmake)
        install(FILES ${MAN1_PATH}/${EXENAME}.1 DESTINATION 
            ${MAN_INSTALL_DIR}/man1)
    endif()
endfunction ()
