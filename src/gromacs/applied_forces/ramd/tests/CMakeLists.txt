gmx_add_unit_test(RAMDAppliedForcesUnitTest ramd_applied_forces-test
    CPP_SOURCE_FILES
        ramd.cpp
        )

target_link_libraries(ramd_applied_forces-test PRIVATE math)
