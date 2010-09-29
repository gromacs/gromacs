IF (GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})
    add_executable(selection-test
                   selectioncollection.cpp selectionoption.cpp
                   test_main.cpp)
    target_link_libraries(selection-test gromacs gmx ${GTEST_LIBRARIES})
    add_test(SelectionUnitTests selection-test)
ENDIF (GTEST_FOUND)
