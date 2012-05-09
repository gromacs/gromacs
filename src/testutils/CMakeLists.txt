include_directories(${GTEST_INCLUDE_DIRS})
set(COMMON_SOURCES datapath.cpp
                   refdata.cpp
                   stringtest.cpp
                   testoptions.cpp)

add_library(testutils STATIC ${COMMON_SOURCES})
set(TESTUTILS_LIBS testutils ${GTEST_LIBRARIES} ${LIBXML2_LIBRARIES})
target_link_libraries(testutils libgromacs ${GTEST_LIBRARIES} ${LIBXML2_LIBRARIES})

set(TESTUTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TESTUTILS_DIR ${TESTUTILS_DIR} PARENT_SCOPE)
set(TESTUTILS_LIBS ${TESTUTILS_LIBS} PARENT_SCOPE)

add_subdirectory(tests)
