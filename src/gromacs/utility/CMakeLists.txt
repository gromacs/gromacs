file(GLOB UTILITY_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${UTILITY_SOURCES} PARENT_SCOPE)

set(UTILITY_PUBLIC_HEADERS
    arrayref.h
    common.h
    errorcodes.h
    exceptions.h
    flags.h
    format.h
    gmxassert.h
    uniqueptr.h)
install(FILES ${UTILITY_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/utility
        COMPONENT development)
