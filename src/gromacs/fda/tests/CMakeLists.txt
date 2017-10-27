set(testname "FDATests")
set(exename "fda-tests")

gmx_add_gtest_executable(
    ${exename}
    LogicallyErrorComparerTest.cpp
    FDATest.cpp
    PairwiseForcesTest.cpp
)

gmx_register_gtest_test(
    ${testname}
    ${exename}
)
