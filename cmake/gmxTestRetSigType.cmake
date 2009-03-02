# - Define macro to check return type of signals (int/void)
#
#  GMX_TEST_RETSIGTYPE(VARIABLE)
#
#  VARIABLE will be set to the return type of signals - "int" or "void"
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_RETSIGTYPE VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for return type of signals")

	# First check without any special flags
        TRY_COMPILE(RETSIGTYPE_INT_OK "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestRetSigType.c")

        if(RETSIGTYPE_INT_OK)
	    MESSAGE(STATUS "Checking for return type of signals - int")			
            set(${VARIABLE} "int" CACHE INTERNAL "Result of test for signal return type" FORCE)
        else(RETSIGTYPE_INT_OK)
            MESSAGE(STATUS "Checking for return type of signals - void")
      	    set(${VARIABLE} "void" CACHE INTERNAL "Result of test for signal return type" FORCE)
      	endif(RETSIGTYPE_INT_OK)
        
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_RETSIGTYPE VARIABLE)



