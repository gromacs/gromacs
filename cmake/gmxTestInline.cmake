# - Define macro to check inline keyword
#
#  GMX_TEST_INLINE(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE VARIABLE)
    IF(NOT DEFINED TEST_${VARIABLE})

        MESSAGE(STATUS "Checking for inline keyword")

	FOREACH(KEYWORD "__inline__" "__inline" "inline")
            IF(NOT TEST_${VARIABLE})
                TRY_COMPILE(TEST_${VARIABLE} "${CMAKE_BINARY_DIR}"    
                            "${CMAKE_SOURCE_DIR}/cmake/TestInline.c"
                            COMPILE_DEFINITIONS "-DTESTINLINEDEF=${KEYWORD}" )
                SET(LAST_INLINE_KEYWORD ${KEYWORD})
	    ENDIF(NOT TEST_${VARIABLE})
        ENDFOREACH(KEYWORD)

        IF(TEST_${VARIABLE})
            SET(${VARIABLE} ${LAST_INLINE_KEYWORD} CACHE INTERNAL "Inline keyword" FORCE)
            MESSAGE(STATUS "Checking for inline keyword - ${LAST_INLINE_KEYWORD}")
        ELSE(TEST_${VARIABLE})
	    SET(${VARIABLE} " " CACHE INTERNAL "Inline keyword" FORCE)
            MESSAGE(STATUS "Checking for inline keyword - not found")
        ENDIF(TEST_${VARIABLE})

    ENDIF(NOT DEFINED TEST_${VARIABLE})        
ENDMACRO(GMX_TEST_INLINE VARIABLE)




