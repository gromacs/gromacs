gmx_add_unit_test(TrajectoryAnalysisUnitTests trajectoryanalysis-test
                  moduletest.cpp
                  angle.cpp
                  select.cpp)

add_executable(test_selection test_selection.cpp)
target_link_libraries(test_selection libgromacs)
