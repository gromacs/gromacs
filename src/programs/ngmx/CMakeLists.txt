if(X11_FOUND)

set(NGMX_COMMON_SOURCE
    xutil.cpp         xmb.cpp           x11.cpp
    xdlghi.cpp        fgrid.cpp         xdlgitem.cpp
    popup.cpp         pulldown.cpp      xdlg.cpp
    manager.cpp       buttons.cpp       nmol.cpp
    nleg.cpp          dialogs.cpp       logo.cpp
    filter.cpp        molps.cpp)

include_directories(${X11_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List of programs with single corresponding *.cpp source file,
# used to create build rules automatically.
#
set(NGMX_PROGRAMS
    ngmx )

foreach(PROG ${NGMX_PROGRAMS})
        add_executable(${PROG} ${PROG}.cpp ${NGMX_COMMON_SOURCE})
        target_link_libraries(${PROG} libgromacs ${GMX_EXTRA_LIBRARIES}
${X11_LIBRARIES})
#${X11_X11_LIB})
        gmx_add_man_page(${PROG})
        set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG)

install(TARGETS ${NGMX_PROGRAMS}
        COMPONENT ngmx
        RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endif(X11_FOUND)
