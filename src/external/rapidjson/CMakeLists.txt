SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

PROJECT(RapidJSON CXX)

option(RAPIDJSON_HAS_STDSTRING "" OFF)
if(RAPIDJSON_HAS_STDSTRING)
    add_definitions(-DRAPIDJSON_HAS_STDSTRING)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_custom_target(travis_doc)
add_custom_command(TARGET travis_doc
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/travis-doxygen.sh)

install(DIRECTORY include/rapidjson
    DESTINATION "${INCLUDE_INSTALL_DIR}"
    COMPONENT dev)

# Provide config and version files to be used by other applications
# ===============================

export(PACKAGE ${PROJECT_NAME})

# cmake-modules
CONFIGURE_FILE(${PROJECT_NAME}Config.cmake.in
    ${PROJECT_NAME}Config.cmake
    @ONLY)
CONFIGURE_FILE(${PROJECT_NAME}ConfigVersion.cmake.in
    ${PROJECT_NAME}ConfigVersion.cmake
    @ONLY)
INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION "${CMAKE_INSTALL_DIR}"
    COMPONENT dev)
