file(GLOB COMMANDLINE_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${COMMANDLINE_SOURCES} PARENT_SCOPE)

set(COMMANDLINE_PUBLIC_HEADERS
    cmdlinehelpwriter.h
    cmdlineparser.h)
install(FILES ${COMMANDLINE_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/commandline
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
