# - Define macro to check if isfinite or _isfinite exists
#
#  GMX_TEST_ISFINITE(VARIABLE)
#
#  VARIABLE will be set to true if ISFINITE exists
#
#  GMX_TEST__ISFINITE(VARIABLE)
#
#  VARIABLE will be set to true if _ISFINITE exists
#

MACRO(GMX_TEST_ISFINITE VARIABLE)
    MESSAGE(STATUS "Checking for ISFINITE")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    set(CMAKE_REQUIRED_LIBRARIES "m")
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  isfinite(f);
}" ISFINITE_COMPILE_OK)

    if(ISFINITE_COMPILE_OK)
        MESSAGE(STATUS "Checking for ISFINITE - yes")
            set(${VARIABLE} ${ISFINITE_COMPILE_OK}
                "Result of test for ISFINITE")
    else(ISFINITE_COMPILE_OK)
        MESSAGE(STATUS "Checking for ISFINITE - no")
    endif(ISFINITE_COMPILE_OK)
ENDMACRO(GMX_TEST_ISFINITE VARIABLE)

MACRO(GMX_TEST__ISFINITE VARIABLE)
    MESSAGE(STATUS "Checking for _ISFINITE")

    set(CMAKE_REQUIRED_INCLUDES "math.h")
    set(CMAKE_REQUIRED_LIBRARIES "m")
    check_c_source_compiles(
      "#include <math.h>
int main(void) {
  float f;
  _isfinite(f);
}" _ISFINITE_COMPILE_OK)

    if(_ISFINITE_COMPILE_OK)
        MESSAGE(STATUS "Checking for _ISFINITE - yes")
            set(${VARIABLE} ${_ISFINITE_COMPILE_OK}
                "Result of test for _ISFINITE")
    else(_ISFINITE_COMPILE_OK)
        MESSAGE(STATUS "Checking for _ISFINITE - no")
    endif(_ISFINITE_COMPILE_OK)
ENDMACRO(GMX_TEST__ISFINITE VARIABLE)
