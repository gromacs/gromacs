add_gtest_or_gmock_test(OptionsUnitTests options-test
                        cmdlineparser.cpp option.cpp optionsassigner.cpp
                        timeunitmanager.cpp
                        GMOCK_ONLY abstractoptionstorage.cpp)
