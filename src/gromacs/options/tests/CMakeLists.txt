IF (GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIR})
    add_executable(options-test
                   cmdlineparser.cpp option.cpp optionsassigner.cpp
                   test_main.cpp)
    target_link_libraries(options-test libgromacs ${GTEST_LIBRARIES})
    add_test(OptionsUnitTests options-test)
ENDIF (GTEST_FOUND)
