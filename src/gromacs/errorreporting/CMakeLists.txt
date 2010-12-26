file(GLOB ERRORREPORTING_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${ERRORREPORTING_SOURCES} PARENT_SCOPE)

set(ERRORREPORTING_PUBLIC_HEADERS
    abstracterrorreporter.h
    emptyerrorreporter.h
    errorcontext.h
    standarderrorreporter.h)
install(FILES ${ERRORREPORTING_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/errorreporting
        COMPONENT development)
