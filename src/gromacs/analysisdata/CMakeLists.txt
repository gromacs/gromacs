file(GLOB ANALYSISDATA_SOURCES *.cpp modules/*.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${ANALYSISDATA_SOURCES} PARENT_SCOPE)

set(ANALYSISDATA_PUBLIC_HEADERS
    abstractdata.h
    analysisdata.h
    arraydata.h
    dataframe.h
    datamodule.h)
install(FILES ${ANALYSISDATA_PUBLIC_HEADERS}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/analysisdata
        COMPONENT development)

add_subdirectory(modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
