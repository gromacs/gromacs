
# - Check the username performing the build, as well as date and time
#
# gmx_set_build_information()
#
# The macro variables will be set to the user/host/cpu used for configuration,
# or anonymous/unknown if it cannot be detected (windows)
#
# BUILD_TIME
# BUILD_USER
# BUILD_HOST
# BUILD_CPU_VENDOR
# BUILD_CPU_BRAND
# BUILD_CPU_FAMILY
# BUILD_CPU_MODEL
# BUILD_CPU_STEPPING
# BUILD_CPU_FEATURES
#

# we rely on inline asm support for GNU!
include(gmxTestInlineASM)

macro(gmx_set_build_information)
    IF(NOT DEFINED BUILD_USER)

    gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

    if(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "-DGMX_X86_GCC_INLINE_ASM")
    else(GMX_X86_GCC_INLINE_ASM)
        set(GCC_INLINE_ASM_DEFINE "")
    endif(GMX_X86_GCC_INLINE_ASM)

    message(STATUS "Setting build user/date/host/cpu information")
    if(CMAKE_HOST_UNIX)
        execute_process( COMMAND date     OUTPUT_VARIABLE TMP_TIME    OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND whoami   OUTPUT_VARIABLE TMP_USER       OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process( COMMAND hostname OUTPUT_VARIABLE TMP_HOSTNAME   OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(BUILD_USER    "@TMP_USER@\@@TMP_HOSTNAME@ [CMAKE]" CACHE INTERNAL "Build user")
        set(BUILD_TIME    "@TMP_TIME@" CACHE INTERNAL "Build date & time")
        execute_process( COMMAND uname -srm OUTPUT_VARIABLE TMP_HOST OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(BUILD_HOST    "@TMP_HOST@" CACHE INTERNAL "Build host & architecture")
        message(STATUS "Setting build user & time - OK")
    else(CMAKE_HOST_UNIX)
        set(BUILD_USER    "Anonymous@unknown [CMAKE]" CACHE INTERNAL "Build user")
        set(BUILD_TIME    "Unknown date" CACHE INTERNAL "Build date & time")
        set(BUILD_HOST    "@CMAKE_HOST_SYSTEM@ @CMAKE_HOST_SYSTEM_PROCESSOR@" CACHE INTERNAL "Build host & architecture")
        message(STATUS "Setting build user & time - not on Unix, using anonymous")
    endif(CMAKE_HOST_UNIX)

    if(NOT CMAKE_CROSSCOMPILING)
        # Get CPU acceleration information
        try_run(GMX_CPUID_RUN_VENDOR GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_VENDOR ARGS "-vendor")
        try_run(GMX_CPUID_RUN_BRAND GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_BRAND ARGS "-brand")
        try_run(GMX_CPUID_RUN_FAMILY GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_FAMILY ARGS "-family")
        try_run(GMX_CPUID_RUN_MODEL GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_MODEL ARGS "-model")
       try_run(GMX_CPUID_RUN_STEPPING GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_STEPPING ARGS "-stepping")
        try_run(GMX_CPUID_RUN_FEATURES GMX_CPUID_COMPILED
            ${CMAKE_BINARY_DIR}
            ${CMAKE_SOURCE_DIR}/src/gmxlib/gmx_cpuid.c
            COMPILE_DEFINITIONS "@GCC_INLINE_ASM_DEFINE@ -I${CMAKE_SOURCE_DIR}/include -DGMX_CPUID_STANDALONE"
            RUN_OUTPUT_VARIABLE OUTPUT_CPU_FEATURES ARGS "-features")

        string(STRIP "@OUTPUT_CPU_VENDOR@" OUTPUT_CPU_VENDOR)
        string(STRIP "@OUTPUT_CPU_BRAND@" OUTPUT_CPU_BRAND)
        string(STRIP "@OUTPUT_CPU_FAMILY@" OUTPUT_CPU_FAMILY)
        string(STRIP "@OUTPUT_CPU_MODEL@" OUTPUT_CPU_MODEL)
        string(STRIP "@OUTPUT_CPU_STEPPING@" OUTPUT_CPU_STEPPING)
        string(STRIP "@OUTPUT_CPU_FEATURES@" OUTPUT_CPU_FEATURES)

        if(GMX_CPUID_RUN_VENDOR EQUAL 0)
            set(BUILD_CPU_VENDOR   "@OUTPUT_CPU_VENDOR@"   CACHE INTERNAL "Build CPU vendor")
        else()
            set(BUILD_CPU_VENDOR   "Unknown, detect failed" CACHE INTERNAL "Build CPU vendor")
        endif()
        if(GMX_CPUID_RUN_BRAND EQUAL 0)
            set(BUILD_CPU_BRAND    "@OUTPUT_CPU_BRAND@"    CACHE INTERNAL "Build CPU brand")
        else()
            set(BUILD_CPU_BRAND    "Unknown, detect failed" CACHE INTERNAL "Build CPU brand")
        endif()
        if(GMX_CPUID_RUN_FAMILY EQUAL 0)
            set(BUILD_CPU_FAMILY   "@OUTPUT_CPU_FAMILY@"   CACHE INTERNAL "Build CPU family")
        else()
            set(BUILD_CPU_FAMILY   "0"                     CACHE INTERNAL "Build CPU family")
        endif()
        if(GMX_CPUID_RUN_MODEL EQUAL 0)
            set(BUILD_CPU_MODEL    "@OUTPUT_CPU_MODEL@"    CACHE INTERNAL "Build CPU model")
        else()
            set(BUILD_CPU_MODEL    "0"                     CACHE INTERNAL "Build CPU model")
        endif()
        if(GMX_CPUID_RUN_STEPPING EQUAL 0)
            set(BUILD_CPU_STEPPING "@OUTPUT_CPU_STEPPING@" CACHE INTERNAL "Build CPU stepping")
        else()
            set(BUILD_CPU_STEPPING "0"                     CACHE INTERNAL "Build CPU stepping")
        endif()
            if(GMX_CPUID_RUN_FEATURES EQUAL 0)
            set(BUILD_CPU_FEATURES "@OUTPUT_CPU_FEATURES@" CACHE INTERNAL "Build CPU features")
        else()
            set(BUILD_CPU_FEATURES ""                      CACHE INTERNAL "Build CPU features")
        endif()

    else(NOT CMAKE_CROSSCOMPILING)
        
        set(BUILD_CPU_VENDOR   "Unknown, cross-compiled"   CACHE INTERNAL "Build CPU vendor")
        set(BUILD_CPU_BRAND    "Unknown, cross-compiled"    CACHE INTERNAL "Build CPU brand")
        set(BUILD_CPU_FAMILY   "0"   CACHE INTERNAL "Build CPU family")
        set(BUILD_CPU_MODEL    "0"    CACHE INTERNAL "Build CPU model")
        set(BUILD_CPU_STEPPING "0" CACHE INTERNAL "Build CPU stepping")
        set(BUILD_CPU_FEATURES "" CACHE INTERNAL "Build CPU features")

    endif(NOT CMAKE_CROSSCOMPILING)

    ENDIF(NOT DEFINED BUILD_USER)
endmacro(gmx_set_build_information)

