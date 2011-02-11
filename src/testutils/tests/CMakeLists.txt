if (TESTUTILS_HAVE_REFDATA)
    add_gtest_test(TestUtilsUnitTests testutils-test
                   refdata.cpp)
endif (TESTUTILS_HAVE_REFDATA)
