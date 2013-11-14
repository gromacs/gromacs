# - Check the username performing the build, as well as date and time
#
# gmx_detect_acceleration(GMX_SUGGESTED_CPU_ACCELERATION)
#
# Try to detect CPU information and suggest an acceleration option
# (such as SSE/AVX) that fits the current CPU. These functions assume
# that gmx_detect_target_architecture() has already been run, so that
# things like GMX_TARGET_X86 are already available.
#
# Sets ${GMX_SUGGESTED_CPU_ACCELERATION} in the parent scope if
# GMX_CPU_ACCELERATION is not set (e.g. by the user, or a previous run
# of CMake).
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

function(gmx_suggest_x86_acceleration _suggested_acceleration)

    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "")
    endif(GMX_X86_GCC_INLINE_ASM)

    message(STATUS "Detecting best acceleration for this CPU")

    # Get CPU acceleration information
    set(_compile_definitions "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/src/gromacs/legacyheaders -DGMX_CPUID_STANDALONE")
    if(GMX_TARGET_X86)
        set(_compile_definitions "${_compile_definitions} -DGMX_TARGET_X86")
    endif()
    try_run(GMX_CPUID_RUN_ACC GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gromacs/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS ${_compile_definitions}
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

    set(${_suggested_acceleration} "@OUTPUT_ACC@" PARENT_SCOPE)
    message(STATUS "Detected best acceleration for this CPU - @OUTPUT_ACC@")
endfunction()

function(gmx_detect_acceleration _suggested_acceleration)
    if(NOT DEFINED GMX_CPU_ACCELERATION)
        if(GMX_TARGET_BGQ)
            set(${_suggested_acceleration} "IBM_QPX")
        elseif(GMX_TARGET_X86)
            gmx_suggest_x86_acceleration(${_suggested_acceleration})
        else()
            set(${_suggested_acceleration} "None")
        endif()

        set(${_suggested_acceleration} ${${_suggested_acceleration}} PARENT_SCOPE)
    endif()
endfunction()
