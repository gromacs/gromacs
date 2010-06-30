
# Test C flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CFLAGSVAR.
MACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_C_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
        IF (${VARIABLE})
            SET (${CFLAGSVAR} "${FLAGS} ${${CFLAGSVAR}}")
        ENDIF (${VARIABLE}) 
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)

# Test C++ flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CXXFLAGSVAR.
MACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_CXX_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
        IF (${VARIABLE})
            SET (${CXXFLAGSVAR} "${FLAGS} ${${CXXFLAGSVAR}}")
        ENDIF (${VARIABLE}) 
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)


# This is the actual exported function to be called 
MACRO(gmx_c_flags)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# check gcc
if(CMAKE_COMPILER_IS_GNUCC)
    GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused" GMXC_CFLAGS)
    # new in gcc 4.5
    GMX_TEST_CFLAG(CFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CFLAGS)
    GMX_TEST_CFLAG(CFLAGS_COPT "-fomit-frame-pointer -finline-functions -funroll-all-loops" GMXC_CFLAGS_RELEASE)
endif()

# check g++
if(CMAKE_COMPILER_IS_GNUCXX )

    GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused" GMXC_CXXFLAGS)
  # new in gcc 4.5
    GMX_TEST_CXXFLAG(CXXFLAGS_EXCESS_PREC "-fexcess-precision=fast" 
                      GMXC_CXXFLAGS)
    GMX_TEST_CXXFLAG(CXXFLAGS_COPT "-fomit-frame-pointer -finline-functions -funroll-all-loops" GMXC_CXXFLAGS_RELEASE)
endif()



# icc
if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if (NOT WIN32) 
        GMX_TEST_CFLAG(CFLAGS_OPT "-ip -w -funroll-all-loops -std=gnu99" 
                        GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_SSE2 "-msse2" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_X86 "-mtune=core2" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_IA64 "-mtune=itanium2" GMXC_CFLAGS_RELEASE)
    else()
        GMX_TEST_CFLAG(CFLAGS_SSE2 "/arch:SSE2" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_X86 "/Qip" GMXC_CFLAGS_RELEASE)
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

    if (NOT WIN32) 
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip -w -funroll-all-loops -std=gnu99" 
                          GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "-msse2" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_X86 "-mtune=core2" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_IA64 "-mtune=itanium2" 
                          GMXC_CXXFLAGS_RELEASE)
    else()
        GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "/arch:SSE2" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_X86 "/Qip" GMXC_CXXFLAGS_RELEASE)
    endif()
endif()

#TODO: other compilers, like xlc

# now actually set the flags:
if (NOT DEFINED GMXCFLAGS_SET)
    set(GMXCFLAGS_SET true CACHE INTERNAL "Whether to reset the C flags" FORCE)
    # C
    set(CMAKE_CFLAGS "${GMXC_CFLAGS} ${CMAKE_CFLAGS}" 
        CACHE STRING "Flags used by the compiler during all build types" FORCE)
    set(CMAKE_CFLAGS_RELEASE "${GMXC_CFLAGS_RELEASE} ${CMAKE_CFLAGS_RELEASE}" 
        CACHE STRING "Flags used by the compiler during release build" FORCE)
    # C++
    set(CMAKE_CXXFLAGS "${GMXC_CXXFLAGS} ${CMAKE_CXXFLAGS}" 
        CACHE STRING "Flags used by the compiler during all build types" FORCE)
    set(CMAKE_CXXFLAGS_RELEASE 
        "${GMXC_CXXFLAGS_RELEASE} ${CMAKE_CXXFLAGS_RELEASE}" 
        CACHE STRING "Flags used by the compiler during all release build" 
        FORCE)
endif()

ENDMACRO(gmx_c_flags)
