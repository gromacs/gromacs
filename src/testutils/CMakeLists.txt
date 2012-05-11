include_directories(${GMOCK_INCLUDE_DIRS})
set(COMMON_SOURCES datapath.cpp
                   mock_helptopic.cpp
                   refdata.cpp
                   stringtest.cpp
                   testoptions.cpp)

add_library(testutils STATIC ${COMMON_SOURCES})
set(TESTUTILS_LIBS testutils ${GMOCK_LIBRARIES} ${LIBXML2_LIBRARIES})
target_link_libraries(testutils libgromacs ${GMOCK_LIBRARIES} ${LIBXML2_LIBRARIES})

set(TESTUTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TESTUTILS_DIR ${TESTUTILS_DIR} PARENT_SCOPE)
set(TESTUTILS_LIBS ${TESTUTILS_LIBS} PARENT_SCOPE)

add_subdirectory(tests)
