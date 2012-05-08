if (TESTUTILS_HAVE_REFDATA)
    add_gtest_test(CommandLineUnitTests commandline-test
                   cmdlinehelpwriter.cpp
                   cmdlineparser.cpp
                   cmdlinetest.cpp)
endif ()
