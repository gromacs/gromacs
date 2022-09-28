if ("$ENV{RUNNER_OS}" STREQUAL "Windows" AND NOT "x$ENV{ENVIRONMENT_SCRIPT}" STREQUAL "x")
  execute_process(
    COMMAND "$ENV{ENVIRONMENT_SCRIPT}" && set
    OUTPUT_FILE environment_script_output.txt
  )
  file(STRINGS environment_script_output.txt output_lines)
  foreach(line IN LISTS output_lines)
    if (line MATCHES "^([a-zA-Z0-9_-]+)=(.*)$")
      set(ENV{${CMAKE_MATCH_1}} "${CMAKE_MATCH_2}")
    endif()
  endforeach()
endif()

set(path_separator ":")
if ("$ENV{RUNNER_OS}" STREQUAL "Windows")
  set(path_separator ";")
endif()
set(ENV{PATH} "$ENV{GITHUB_WORKSPACE}${path_separator}$ENV{PATH}")

message(STATUS "Using GPU_VAR: $ENV{GPU_VAR}")

execute_process(
  COMMAND cmake
    -S .
    -B build
    -D CMAKE_BUILD_TYPE=$ENV{BUILD_TYPE}
    -G Ninja
    -D CMAKE_MAKE_PROGRAM=ninja
    -D CMAKE_C_COMPILER_LAUNCHER=ccache
    -D CMAKE_CXX_COMPILER_LAUNCHER=ccache
    -D GMX_COMPILER_WARNINGS=ON
    -D GMX_DEFAULT_SUFFIX=OFF
    -D GMX_GPU=$ENV{GPU_VAR}
    -D GMX_SIMD=None
    -D GMX_FFT_LIBRARY=FFTPACK
    -D GMX_OPENMP=$ENV{OPENMP_VAR}
    -D REGRESSIONTEST_DOWNLOAD=ON
  RESULT_VARIABLE result
)
if (NOT result EQUAL 0)
  message(FATAL_ERROR "Bad exit status")
endif()
