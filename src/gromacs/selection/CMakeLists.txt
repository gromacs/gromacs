file(GLOB SELECTION_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${SELECTION_SOURCES} PARENT_SCOPE)

set(SELECTION_PUBLIC_HEADERS
    centerofmass.h
    nbsearch.h
    poscalc.h
    indexutil.h
    position.h
    selection.h
    selectionenums.h
    selectionoption.h
    selparam.h
    selmethod.h
    selvalue.h)
install(FILES ${SELECTION_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/selection
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
