execute_process(
  COMMAND ctest --output-on-failure
  WORKING_DIRECTORY build
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE output
  ECHO_OUTPUT_VARIABLE ECHO_ERROR_VARIABLE
)
if (NOT result EQUAL 0)
  string(REGEX MATCH "[0-9]+% tests.*[0-9.]+ sec.*$" test_results "${output}")
  string(REPLACE "\n" "%0A" test_results "${test_results}")
  message("::error::${test_results}")
  message(FATAL_ERROR "Running tests failed!")
endif()
