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



# - Define macro to check if SIGUSR1 is defined
#
#  GMX_TEST_SIGUSR1(VARIABLE)
#
#  VARIABLE will be set if SIGUSR1 is present in signal.h
#

MACRO(GMX_TEST_SIGUSR1 VARIABLE)
    IF(NOT DEFINED HAVE_${VARIABLE})
        
        MESSAGE(STATUS "Checking for SIGUSR1")

        TRY_COMPILE(HAVE_${VARIABLE} "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestSIGUSR1.c")

	IF(HAVE_${VARIABLE})	    
            MESSAGE(STATUS "Checking for SIGUSR1 - found")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for SIGUSR1" FORCE)
        ELSE(HAVE_${VARIABLE})
            MESSAGE(STATUS "Checking for SIGUSR1 - not found")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for SIGUSR1" FORCE)
        ENDIF(HAVE_${VARIABLE})
        
    ENDIF(NOT DEFINED HAVE_${VARIABLE})
ENDMACRO(GMX_TEST_SIGUSR1 VARIABLE)


