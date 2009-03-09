# - Define macro to check if pipes are supported
#
#  GMX_TEST_PIPES(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_PIPES VARIABLE)
    IF(NOT DEFINED ${VARIABLE})
        
        MESSAGE(STATUS "Checking for pipe support")

        TRY_COMPILE(HAVE_PIPES "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestPipes.c")

    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_PIPES VARIABLE)




