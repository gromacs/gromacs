file(GLOB UTILITY_SOURCES *.cpp)
list(APPEND LIBGROMACS_SOURCES ${UTILITY_SOURCES})

set(GENERATED_HEADER_CONFIG ${CMAKE_CURRENT_BINARY_DIR}/gmx_header_config_gen.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gmx_header_config_gen.h.cmakein
               ${GENERATED_HEADER_CONFIG})
# Add target that generates gitversion.c every time make is run
# if git version info is requested
if (GMX_GIT_VERSION_INFO)
    set(GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/gitversion.c)
    add_custom_target(gmx_version ALL
            COMMAND ${CMAKE_COMMAND}
                -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                -D GIT_VERSION="${GIT_VERSION}"
                -D PROJECT_VERSION="${PROJECT_VERSION}"
                -D PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
                -D VERSION_C_CMAKEIN="${CMAKE_CURRENT_SOURCE_DIR}/gitversion.c.cmakein"
                -D VERSION_C_OUT=${GENERATED_VERSION_FILE}
                -P ${CMAKE_SOURCE_DIR}/cmake/gmxGenerateVersionInfo.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gitversion.c.cmakein
            COMMENT "Generating git version information")
    set_source_files_properties(${GENERATED_VERSION_FILE}
                                PROPERTIES GENERATED true)
    list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE})
endif()
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} PARENT_SCOPE)

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
install(FILES ${GENERATED_HEADER_CONFIG}
        DESTINATION ${INCL_INSTALL_DIR}/gromacs/utility
        COMPONENT development)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
