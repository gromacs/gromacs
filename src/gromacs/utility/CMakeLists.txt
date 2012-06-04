file(GLOB UTILITY_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${UTILITY_SOURCES} PARENT_SCOPE)

# Needs to be created in the source tree to allow relative paths in headers to
# be resolved correctly before installation
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gmx_header_config.h.cmakein
               ${CMAKE_CURRENT_SOURCE_DIR}/gmx_header_config.h)

set(UTILITY_PUBLIC_HEADERS
    arrayref.h
    common.h
    errorcodes.h
    exceptions.h
    file.h
    flags.h
    gmx_header_config.h
    gmxassert.h
    programinfo.h
    stringutil.h
    uniqueptr.h)
install(FILES ${UTILITY_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/utility
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
