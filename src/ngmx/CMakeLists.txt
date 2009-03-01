
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


add_executable(ngmx ngmx.c ${NGMX_COMMON_SOURCE})
target_link_libraries(ngmx md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(xrama xrama.c ${NGMX_COMMON_SOURCE})
target_link_libraries(xrama md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(highway highway.c ${NGMX_COMMON_SOURCE})
target_link_libraries(highway md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(gmxlogo gmxlogo.c ${NGMX_COMMON_SOURCE})
target_link_libraries(gmxlogo md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(dlg dlg.c ${NGMX_COMMON_SOURCE})
target_link_libraries(dlg md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(scrollw scrollw.c ${NGMX_COMMON_SOURCE})
target_link_libraries(scrollw md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

add_executable(showcol showcol.c ${NGMX_COMMON_SOURCE})
target_link_libraries(showcol md ${X11_LIBRARIES} ${GROMACS_EXTRA_LIBRARIES})

endif(X11_FOUND)
