if (TESTUTILS_HAVE_REFDATA)
    add_gtest_or_gmock_test(CommandLineUnitTests commandline-test
                            cmdlinehelpwriter.cpp
                            cmdlineparser.cpp
                            cmdlinetest.cpp
                            GMOCK_ONLY
                            cmdlinemodulemanager.cpp)
endif ()
