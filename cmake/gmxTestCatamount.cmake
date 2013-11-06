# - Define function to check we are compiling for CRAY XT catamount
#
#  GMX_TEST_CATAMOUNT(VARIABLE)
#
#  VARIABLE will be set to true if we are compiling for catamount
#

function(GMX_TEST_CATAMOUNT VARIABLE)
    include(CheckCSourceCompiles)
    check_c_source_compiles(
        "int
main()
{
#if defined __QK_USER__
  return 0;
#else
#  error not catamount
#endif
}" CATAMOUNT_COMPILE_OK)
    set(${VARIABLE} ${CATAMOUNT_COMPILE_OK} PARENT_SCOPE)
endfunction()
