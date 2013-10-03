# - Define macro to check we are compiling for CRAY XT catamount
#
#  GMX_TEST_CATAMOUNT(VARIABLE)
#
#  VARIABLE will be set to true if we are compiling for catamount
#

MACRO(GMX_TEST_CATAMOUNT VARIABLE)
    TRY_COMPILE(CATAMOUNT_COMPILE_OK "${CMAKE_BINARY_DIR}"
        "${CMAKE_SOURCE_DIR}/cmake/TestCatamount.c")
    set(${VARIABLE} ${CATAMOUNT_COMPILE_OK} PARENT_SCOPE)
ENDMACRO(GMX_TEST_CATAMOUNT VARIABLE)
