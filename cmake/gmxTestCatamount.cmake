# - Define macro to check we are compiling for CRAY XT catamount
#
#  GMX_TEST_CATAMOUNT(VARIABLE)
#
#  VARIABLE will be set to true if we are compiling for catamount
#

MACRO(GMX_TEST_CATAMOUNT VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for CRAY XT Catamount compile")

	# First check without any special flags
        TRY_COMPILE(CATAMOUNT_COMPILE_OK "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestCatamount.c")

        if(CATAMOUNT_COMPILE_OK)
	    MESSAGE(STATUS "Checking for CRAY XT Catamount target - yes")			
        else(CATAMOUNT_COMPILE_OK)
            MESSAGE(STATUS "Checking for CRAY XT Catamount target - no")
      	endif(CATAMOUNT_COMPILE_OK)

        set(${VARIABLE} ${CATAMOUNT_COMPILE_OK} CACHE INTERNAL 
            "Result of test for CRAY XT Catamount target" FORCE)
        
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_CATAMOUNT VARIABLE)



