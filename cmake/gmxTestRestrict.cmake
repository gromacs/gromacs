# - Define macro to check restrict keyword
#
#  GMX_TEST_RESTRICT(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_RESTRICT VARIABLE)
    IF(NOT DEFINED TEST_${VARIABLE})

        MESSAGE(STATUS "Checking for restrict keyword")

# Start with __restrict__, since that is the C++ default keyword.
	FOREACH(KEYWORD "__restrict__" "__restrict" "restrict")
            IF(NOT TEST_${VARIABLE})
                TRY_COMPILE(TEST_${VARIABLE} "${CMAKE_BINARY_DIR}"    
                            "${CMAKE_SOURCE_DIR}/cmake/TestRestrict.c"
                            COMPILE_DEFINITIONS "-DTESTRESTRICTDEF=${KEYWORD}" )
                SET(CHK_RESTRICT_KEYWORD ${KEYWORD})
            ENDIF(NOT TEST_${VARIABLE})
        ENDFOREACH(KEYWORD)
             
        IF(TEST_${VARIABLE})
            SET(${VARIABLE} ${KEYWORD})
            MESSAGE(STATUS "Checking for restrict keyword - ${CHK_RESTRICT_KEYWORD}")
        ELSE(TEST_${VARIABLE})
	    SET(${VARIABLE} " ")
            MESSAGE(STATUS "Checking for restrict keyword - not found")
        ENDIF(TEST_${VARIABLE})

    ENDIF(NOT DEFINED TEST_${VARIABLE})        
ENDMACRO(GMX_TEST_RESTRICT VARIABLE)




