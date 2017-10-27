set(testname "regressiontests/fda-analysis")
set(exename "regressiontests-fda-analysis")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests 
    FDAGetStressTest.cpp
    FDAGraphTest.cpp
    FDAShortestPathTest.cpp
    FDAViewStressTest.cpp
    PDBTest.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
)

target_compile_definitions(${exename} PRIVATE REGRESSIONTEST_PATH=${REGRESSIONTEST_PATH})

gmx_register_gtest_test(
    ${testname}
    ${exename}
    INTEGRATION_TEST
)
