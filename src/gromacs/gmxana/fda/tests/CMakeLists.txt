set(testname "FDAGraphTest")
set(exename "fda-graph-test")

gmx_add_gtest_executable(
    ${exename}
    BoostGraphTest.cpp
)

gmx_register_gtest_test(
    ${testname}
    ${exename}
)
