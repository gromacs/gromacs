PROJECT(RapidJSON CXX)

option(RAPIDJSON_HAS_STDSTRING "" OFF)
if(RAPIDJSON_HAS_STDSTRING)
    add_definitions(-DRAPIDJSON_HAS_STDSTRING)
endif()

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#add_custom_target(travis_doc)
#add_custom_command(TARGET travis_doc
#    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/travis-doxygen.sh)

install(DIRECTORY include/rapidjson
    DESTINATION "${INCLUDE_INSTALL_DIR}"
    COMPONENT dev)
#    )


