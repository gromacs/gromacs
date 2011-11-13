if (TESTUTILS_HAVE_REFDATA)
    add_gmock_test(AnalysisDataUnitTests analysisdata-test
                   analysisdata.cpp arraydata.cpp
                   average.cpp datatest.cpp histogram.cpp
                   mock_module.cpp)
endif ()
