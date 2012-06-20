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
    try_run(GMX_DETECTCPU_RUN_ACC GMX_DETECTCPU_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_detectcpu.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_DETECTCPU_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_TMP ARGS "-acceleration")
    string(STRIP "@OUTPUT_TMP@" OUTPUT_ACC)

    message(STATUS "Detecting best acceleration for this CPU - @OUTPUT_ACC@")

    set(${GMX_SUGGESTED_ACCELERATION}    "@OUTPUT_ACC@" CACHE INTERNAL "Gromacs CPU Acceleration")

    ENDIF(NOT DEFINED ${GMX_SUGGESTED_ACCELERATION})
endmacro(gmx_detect_acceleration GMX_SUGGESTED_ACCELERATION)

