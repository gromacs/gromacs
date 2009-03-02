# - Define macro to check if system XDR routines are present
#
#  GMX_TEST_XDR(VARIABLE)
#
#  VARIABLE will be set to true if XDR support is present
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_XDR VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for system XDR support")

	# First check without any special flags
        TRY_COMPILE(XDR_COMPILE_OK "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestXDR.c")

        if(XDR_COMPILE_OK)
	    MESSAGE(STATUS "Checking for system XDR support - present")			
        else(XDR_COMPILE_OK)
            MESSAGE(STATUS "Checking for system XDR support - not present")
      	endif(XDR_COMPILE_OK)

        set(${VARIABLE} ${XDR_COMPILE_OK} CACHE INTERNAL "Result of test for system XDR support" FORCE)
        
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_XDR VARIABLE)



