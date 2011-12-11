file(GLOB TRAJECTORYANALYSIS_SOURCES *.cpp modules/*.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${TRAJECTORYANALYSIS_SOURCES} PARENT_SCOPE)

set(TRAJECTORYANALYSIS_PUBLIC_HEADERS
    analysismodule.h
    analysissettings.h
    cmdlinerunner.h)
install(FILES ${TRAJECTORYANALYSIS_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/trajectoryanalysis
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
