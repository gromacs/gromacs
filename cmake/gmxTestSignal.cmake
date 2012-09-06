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


