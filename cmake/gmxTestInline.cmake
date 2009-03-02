# - Define macro to check inline keyword
#
#  GMX_TEST_INLINE(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE VARIABLE)
    IF(NOT DEFINED ${VARIABLE})
        
        MESSAGE(STATUS "Checking for inline keyword")

	FOREACH(KEYWORD "inline" "__inline__" "__inline")
            IF(NOT DEFINED C_INLINE)
                TRY_COMPILE(C_HAS_${KEYWORD} "${CMAKE_BINARY_DIR}"    
                            "${CMAKE_SOURCE_DIR}/cmake/TestInline.c")
                IF(C_HAS_${KEYWORD})
                    SET(C_INLINE TRUE)
                    SET(${VARIABLE} ${KEYWORD})
                    MESSAGE(STATUS "Checking for inline keyword - ${KEYWORD}")
                ENDIF(C_HAS_${KEYWORD})
            ENDIF(NOT DEFINED C_INLINE)
        ENDFOREACH(KEYWORD)
     
        IF(NOT DEFINED C_INLINE)
            SET(${VARIABLE} " ")
            MESSAGE(STATUS "Checking for inline keyword - not found")
        ENDIF(NOT DEFINED C_INLINE)
        
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_INLINE VARIABLE)




