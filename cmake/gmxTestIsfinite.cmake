# - Define macro to check if isfinite or _isfinite exists
#
#  gmx_test_isfinite(VARIABLE)
#
#  VARIABLE will be set to true if isfinite exists
#
#  gmx_test__isfinite(VARIABLE)
#
#  VARIABLE will be set to true if _isfinite exists
#
#  gmx_test__finite(VARIABLE) - disabled since it doesn't seem to work the way the MSVC docs suggest
#
#  VARIABLE will be set to true if _finite exists
#

MACRO(gmx_test_isfinite VARIABLE)
  if(NOT DEFINED isfinite_compile_ok)
    MESSAGE(STATUS "Checking for isfinite")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    if (HAVE_LIBM)
        set(CMAKE_REQUIRED_LIBRARIES "m")
    endif()
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  isfinite(f);
}" isfinite_compile_ok)

    if(isfinite_compile_ok)
        MESSAGE(STATUS "Checking for isfinite - yes")
    else(isfinite_compile_ok)
        MESSAGE(STATUS "Checking for isfinite - no")
    endif(isfinite_compile_ok)
    set(isfinite_compile_ok "${isfinite_compile_ok}" CACHE INTERNAL "Result of isfinite check")
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(NOT DEFINED isfinite_compile_ok)
  if(isfinite_compile_ok)
    set(${VARIABLE} ${isfinite_compile_ok}
                "Result of test for isfinite")
  endif()
ENDMACRO(gmx_test_isfinite VARIABLE)

MACRO(gmx_test__isfinite VARIABLE)
  if(NOT DEFINED _isfinite_compile_ok)
    MESSAGE(STATUS "Checking for _isfinite")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    if (HAVE_LIBM)
        set(CMAKE_REQUIRED_LIBRARIES "m")
    endif()
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  _isfinite(f);
}" _isfinite_compile_ok)

    if(_isfinite_compile_ok)
        MESSAGE(STATUS "Checking for _isfinite - yes")
    else(_isfinite_compile_ok)
        MESSAGE(STATUS "Checking for _isfinite - no")
    endif(_isfinite_compile_ok)
    set(_isfinite_compile_ok "${_isfinite_compile_ok}" CACHE INTERNAL "Result of _isfinite check")
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(NOT DEFINED _isfinite_compile_ok)
  if(_isfinite_compile_ok)
    set(${VARIABLE} ${_isfinite_compile_ok}
                "Result of test for _isfinite")
  endif()
ENDMACRO(gmx_test__isfinite VARIABLE)

# Necessary for MSVC
MACRO(gmx_test__finite VARIABLE)
  if(NOT DEFINED _finite_compile_ok)
    MESSAGE(STATUS "Checking for _finite")

    set(CMAKE_REQUIRED_INCLUDES "float.h")
    check_c_source_compiles(
      "#include <float.h>
int main(void) {
  float f;
  _finite(f);
}" _finite_compile_ok)

    if(_finite_compile_ok)
        MESSAGE(STATUS "Checking for _finite - yes")
    else(_finite_compile_ok)
        MESSAGE(STATUS "Checking for _finite - no")
    endif(_finite_compile_ok)
    set(_finite_compile_ok "${_finite_compile_ok}" CACHE INTERNAL "Result of _finite check")
    set(CMAKE_REQUIRED_INCLUDES)
  endif(NOT DEFINED _finite_compile_ok)
  if(_finite_compile_ok)
    set(${VARIABLE} ${_finite_compile_ok}
                "Result of test for _finite")
  endif()
ENDMACRO(gmx_test__finite VARIABLE)
