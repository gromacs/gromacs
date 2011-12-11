if (TESTUTILS_HAVE_REFDATA)
    add_gtest_test(SelectionUnitTests selection-test
                   selectioncollection.cpp selectionoption.cpp)
endif ()
