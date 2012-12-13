file(GLOB QMMM_SOURCES *.cpp modules/*.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${QMMM_SOURCES} PARENT_SCOPE)

set(QMMM_PUBLIC_HEADERS
    qmmm.h)
install(FILES ${QMMM_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/qmmm
        COMPONENT development)

#if (BUILD_TESTING)
#    add_subdirectory(tests)
#endif (BUILD_TESTING)
