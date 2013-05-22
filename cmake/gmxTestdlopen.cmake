# - Define macro to check if DLOPEN is defined
#
#  GMX_TEST_DLOPEN(VARIABLE)
#
#  VARIABLE will be set if dlopen is present in dlfcn.h
#

MACRO(GMX_TEST_DLOPEN VARIABLE)
  IF(NOT DEFINED ${VARIABLE})
    MESSAGE(STATUS "Checking for dlopen")

    set(CMAKE_REQUIRED_INCLUDES "dlfcn.h")
    set(CMAKE_REQUIRED_LIBRARIES "dl")
    check_c_source_compiles(
      "#include <dlfcn.h>
int main(void) {
  dlopen(0,0);
}" ${VARIABLE})

    IF(${VARIABLE})
      MESSAGE(STATUS "Checking for dlopen - found")
      set(${VARIABLE} 1 CACHE INTERNAL "Result of test for dlopen" FORCE)
    ELSE()
      MESSAGE(STATUS "Checking for dlopen - not found")
      set(${VARIABLE} 0 CACHE INTERNAL "Result of test for dlopen" FORCE)
    ENDIF()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_LIBRARIES)
  ENDIF()
ENDMACRO()