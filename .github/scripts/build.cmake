if ("$ENV{RUNNER_OS}" STREQUAL "Windows" AND NOT "x$ENV{ENVIRONMENT_SCRIPT}" STREQUAL "x")
  file(STRINGS environment_script_output.txt output_lines)
  foreach(line IN LISTS output_lines)
    if (line MATCHES "^([a-zA-Z0-9_-]+)=(.*)$")
      set(ENV{${CMAKE_MATCH_1}} "${CMAKE_MATCH_2}")
    endif()
  endforeach()
endif()

file(TO_CMAKE_PATH "$ENV{GITHUB_WORKSPACE}" ccache_basedir)
set(ENV{CCACHE_BASEDIR} "${ccache_basedir}")
set(ENV{CCACHE_DIR} "${ccache_basedir}/.ccache")
set(ENV{CCACHE_COMPRESS} "true")
set(ENV{CCACHE_COMPRESSLEVEL} "6")
set(ENV{CCACHE_MAXSIZE} "600M")

execute_process(COMMAND ccache -p)
execute_process(COMMAND ccache -z)

execute_process(
  COMMAND cmake --build build
  RESULT_VARIABLE result-build
  OUTPUT_VARIABLE output-build
  ERROR_VARIABLE output-build
  ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
)
execute_process(
  COMMAND cmake --build build --target tests
  RESULT_VARIABLE result-build-test
  OUTPUT_VARIABLE output-build-test
  ERROR_VARIABLE output-build-test
  ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
)
if ((NOT result-build EQUAL 0) OR (NOT result-build-test EQUAL 0))
  string(REGEX MATCH "FAILED:.*$" error_message_build "${output-build}")
  string(REGEX MATCH "FAILED:.*$" error_message_build_test "${output-build-test}")
  string(REPLACE "\n" "%0A" error_message_build "${error_message_build}")
  string(REPLACE "\n" "%0A" error_message_build_test "${error_message_build_test}")
  message("::error::${error_message_build}")
  message("::error::${error_message_build_test}")
  message(FATAL_ERROR "Build failed")
endif()

execute_process(COMMAND ccache -s)