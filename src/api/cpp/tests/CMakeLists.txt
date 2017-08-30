# For an executable target already defined, CMake tests are added with the
# following syntax:
#
# add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
#            [WORKING_DIRECTORY dir]
#            COMMAND <command> [arg1 [arg2 ...]])
#
# The Gromacs convention to get the appropriate compilation environment for
# unit tests is the gmx_add_unit_test macro, defined in
# src/testutils/TestMacros.cmake, which does not allow arbitrary argument
# passing to add_test. Instead, use set_tests_properties().

# Copy test files to test directory
configure_file(data/topol.tpr topol.tpr COPYONLY)

configure_file(test_runner.in.cpp test_runner.cpp)

gmx_add_unit_test(GmxapiExternalInterfaceTests gmxapi-test
                  test_md.cpp
                  test_runner.cpp
                  test_status.cpp
                  test_system.cpp
                  test_version_helpers.cpp
                  )

set_tests_properties(GmxapiExternalInterfaceTests PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_compile_definitions(gmxapi-test PUBLIC
    GMXAPI_MAJOR=${GMXAPI_MAJOR}
    GMXAPI_MINOR=${GMXAPI_MINOR}
    GMXAPI_PATCH=${GMXAPI_PATCH})
target_link_libraries(gmxapi-test Gromacs::gmxapi)

#add_executable(test_apirunner test_basic_mdrunner.cpp)
#gmx_add_unit_test(gmxapi_runner test_apirunner test_basic_mdrunner.cpp)
#target_link_libraries(test_apirunner gmxapi)
#
#gmx_add_unit_test(gmxapi_context test_apicontext test_context.cpp ../context.cpp ../gmxapi.cpp)
#target_include_directories(test_apicontext PRIVATE $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>)
#target_link_libraries(test_apicontext libgromacs)
#
#gmx_add_unit_test(gmxapi_system test_apisystem test_system.cpp ../system.cpp ../gmxapi.cpp ../atoms.cpp ../md.cpp)
#target_include_directories(test_apisystem PRIVATE $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>)
#target_link_libraries(test_apisystem libgromacs)
#
#gmx_add_unit_test(gmxapi_md test_apimd test_md.cpp)
#target_include_directories(test_apimd PRIVATE $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>)
#target_link_libraries(test_apimd gmxapi)
#
#gmx_add_unit_test(gmxapi_runner test_apirunner test_runner.cpp ../runner.cpp ../md.cpp ../gmxapi.cpp $<TARGET_OBJECTS:mdrun_objlib>)
#target_include_directories(test_apirunner PRIVATE $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>)
#target_link_libraries(test_apirunner libgromacs)
#
##gmx_add_unit_test(gmxapi_basictypes test_apibasictypes test_basictypes.cpp)
#
#gmx_add_unit_test(gmxapi_implementation test_impl test_md_impl.cpp ../md.cpp)
#target_include_directories(test_apirunner PRIVATE $<TARGET_PROPERTY:gmxapi,INCLUDE_DIRECTORIES>)
#target_link_libraries(test_impl libgromacs)
