# - Check the username performing the build, as well as date and time
#
# GMX_DETECT_ACCELERATION(GMX_SUGGESTED_ACCELERATION)
#
# Try to detect CPU information and suggest an acceleration option
# (such as SSE/AVX) that fits the current CPU.
#
# GMX_SUGGESTED_ACCELERATION
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

macro(gmx_detect_acceleration GMX_SUGGESTED_ACCELERATION)
    IF(NOT DEFINED ${GMX_SUGGESTED_ACCELERATION})

    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "")
    endif(GMX_X86_GCC_INLINE_ASM)

    message(STATUS "Detecting best acceleration for this CPU")

    # Get CPU acceleration information
    try_run(GMX_CPUID_RUN_ACC GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/src/gromacs/legacyheaders/ -DGMX_DETECTCPU_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_TMP
            COMPILE_OUTPUT_VARIABLE GMX_CPUID_COMPILE_OUTPUT
            ARGS "-acceleration")

    if(NOT GMX_CPUID_COMPILED)
        message(WARNING "Cannot compile CPUID code, which means no CPU-specific acceleration.")
        message(STATUS "Compile output: ${GMX_CPUID_COMPILE_OUTPUT}")
        set(OUTPUT_TMP "None")
    elseif(NOT GMX_CPUID_RUN_ACC EQUAL 0)
        message(WARNING "Cannot run CPUID code, which means no CPU-specific optimization.")
        message(STATUS "Run output: ${OUTPUT_TMP}")
        set(OUTPUT_TMP "None")
    endif(NOT GMX_CPUID_COMPILED)

    string(STRIP "@OUTPUT_TMP@" OUTPUT_ACC)

    message(STATUS "Detecting best acceleration for this CPU - @OUTPUT_ACC@")

    set(${GMX_SUGGESTED_ACCELERATION}    "@OUTPUT_ACC@" CACHE INTERNAL "GROMACS CPU-specific acceleration")

    ENDIF(NOT DEFINED ${GMX_SUGGESTED_ACCELERATION})
endmacro(gmx_detect_acceleration GMX_SUGGESTED_ACCELERATION)

