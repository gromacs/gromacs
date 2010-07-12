# - Define macro to check GCC x86 inline ASM support
#
#  GMX_TEST_INLINE_ASM_GCC_X86(VARIABLE)
#
#  VARIABLE will be set to true if GCC x86 inline asm works.
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)
    IF(NOT DEFINED ${VARIABLE})
        
        MESSAGE(STATUS "Checking for GCC x86 inline asm")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_gcc_x86.c")

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for GCC x86 inline asm - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
        else(${VARIABLE})
            MESSAGE(STATUS "Checking for GCC x86 inline asm - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for GCC x86 inline asm" FORCE)
      	endif(${VARIABLE})


    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_INLINE_ASM_GCC_X86 VARIABLE)



# - Define macro to check MSVC x86 inline ASM support
#
#  GMX_TEST_INLINE_ASM_MSVC_X86(VARIABLE)
#
#  VARIABLE will be set to true if MSVC x86 inline asm works.
#
#  Remember to have a cmakedefine for it too...

MACRO(GMX_TEST_INLINE_ASM_MSVC_X86 VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for MSVC x86 inline asm")

        TRY_COMPILE(${VARIABLE} "${CMAKE_BINARY_DIR}"
                    "${CMAKE_SOURCE_DIR}/cmake/TestInlineASM_msvc_x86.c")

        if(${VARIABLE})
            MESSAGE(STATUS "Checking for MSVC x86 inline asm - supported")
            set(${VARIABLE} 1 CACHE INTERNAL "Result of test for MSVC x86 inline asm" FORCE)
      	else(${VARIABLE})
            MESSAGE(STATUS "Checking for MSVC x86 inline asm - not supported")
            set(${VARIABLE} 0 CACHE INTERNAL "Result of test for MSVC x86 inline asm" FORCE)
        endif(${VARIABLE})

    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_INLINE_ASM_MSVC_X86 VARIABLE)





