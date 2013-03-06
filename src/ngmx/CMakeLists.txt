if(X11_FOUND)

set(NGMX_COMMON_SOURCE
    xutil.c         xmb.c           x11.c
    xdlghi.c        fgrid.c         xdlgitem.c
    popup.c         pulldown.c      xdlg.c
    manager.c       buttons.c       nmol.c
    nleg.c          dialogs.c       logo.c
    filter.c        molps.c)

include_directories(${X11_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List of programs with single corresponding *.c source file,
# used to create build rules automatically.
#
set(NGMX_PROGRAMS 
    ngmx g_xrama)

foreach(PROG ${NGMX_PROGRAMS})
        add_executable(${PROG} ${PROG}.c ${NGMX_COMMON_SOURCE})
        target_link_libraries(${PROG} gmx ${GMX_EXTRA_LIBRARIES} ${X11_LIBRARIES})        
        set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG) 

list(APPEND TARGETS_NEEDED_BY_MANUAL_BUILD ${NGMX_PROGRAMS})
# Set global variable from the local variable's value
set(TARGETS_NEEDED_BY_MANUAL_BUILD ${TARGETS_NEEDED_BY_MANUAL_BUILD} PARENT_SCOPE)

install(TARGETS ${NGMX_PROGRAMS}
        COMPONENT ngmx
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(X11_FOUND)
