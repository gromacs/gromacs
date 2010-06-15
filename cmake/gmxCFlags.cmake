MACRO(gmx_c_flags)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

if(CMAKE_COMPILER_IS_GNUCC)
# this one is set by cmake depending on the build target
#  CHECK_C_COMPILER_FLAG( "-O3" XFLAGS_O3)
#  if (XFLAGS_O3)
#    set(GMXC_C_FLAGS " -O3 ${GMXC_C_FLAGS}")
#  endif()

  CHECK_C_COMPILER_FLAG( "-Wall -Wno-unused" XFLAGS_WARN)
  if (XFLAGS_WARN)
    set(GMXC_C_FLAGS "-Wall -Wno-unused ${GMXC_C_FLAGS}")
  endif()

  CHECK_C_COMPILER_FLAG( "-std=gnu99" XFLAGS_GNU99)
  if (XFLAGS_GNU99)
    set(GMXC_C_FLAGS "-std=gnu99 ${GMXC_C_FLAGS}")
  endif()

  CHECK_C_COMPILER_FLAG( "-march=native" XFLAGS_MARCH)
  if (XFLAGS_MARCH)
    set(GMXC_C_FLAGS "-march=native ${GMXC_C_FLAGS}")
  endif()

  # new in gcc 4.5
  CHECK_C_COMPILER_FLAG( "-fexcess-precision=fast" XFLAGS_EXCESS_PRECISION)
  if (XFLAGS_EXCESS_PRECISION)
    set(GMXC_C_FLAGS "-fexcess-precision=fast ${GMXC_C_FLAGS}")
  endif()

  CHECK_C_COMPILER_FLAG( "-fomit-frame-pointer -finline-functions -funroll-all-loops" XFLAGS_OPT)
  if (XFLAGS_OPT)
    set(GMXC_C_FLAGS_RELEASE "-fomit-frame-pointer -finline-functions -funroll-all-loops ${GMXC_C_FLAGS_RELEASE}")
  endif()
endif()


if(CMAKE_COMPILER_IS_GNUCXX )
# this one is set by cmake depending on the build target
#  CHECK_CXX_COMPILER_FLAG( "-O3" XXFLAGS_O3)
#  if (XXFLAGS_O3)
##    set(GMXC_CXX_FLAGS " -O3 ${GMXC_CXX_FLAGS}")
#  endif()

  CHECK_CXX_COMPILER_FLAG( "-Wall -Wno-unused" XXFLAGS_WARN)
  if (XXFLAGS_WARN)
    set(GMXC_CXX_FLAGS "-Wall -Wno-unused ${GMXC_CXX_FLAGS}")
  endif()

  CHECK_CXX_COMPILER_FLAG( "-march=native" XXFLAGS_MARCH)
  if (XXFLAGS_MARCH)
    set(GMXC_CXX_FLAGS "-march=native ${GMXC_CXX_FLAGS}")
  endif()
  # new in gcc 4.5
  CHECK_CXX_COMPILER_FLAG( "-fexcess-precision=fast" XXFLAGS_EXCESS_PRECISION)
  if (XXFLAGS_EXCESS_PRECISION)
    set(GMXC_CXX_FLAGS "-fexcess-precision=fast ${GMXC_CXX_FLAGS}")
  endif()

  CHECK_CXX_COMPILER_FLAG( "-fomit-frame-pointer -finline-functions -funroll-all-loops" XXFLAGS_OPT)
  if (XXFLAGS_OPT)
    set(GMXC_CXX_FLAGS_RELEASE "-fomit-frame-pointer -finline-functions -funroll-all-loops ${GMXC_CXX_FLAGS_RELEASE}")
  endif()
endif()


# icc
if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if (NOT WIN32) 
        CHECK_C_COMPILER_FLAG("-ip -w -funroll-all-loops -std=gnu99" XFLAGS_OPT)
        if (XFLAGS_OPT)
            set(GMXC_C_FLAGS "-ip -w -funroll-all-loops -std=gnu99 ${GMXC_C_FLAGS}")
        endif()
        CHECK_C_COMPILER_FLAG("-msse2" XFLAGS_OPT_SSE2)
        if (XFLAGS_OPT_SSE2)
            set(GMXC_C_FLAGS "-msse2 ${GMXC_C_FLAGS}")
        endif()
        CHECK_C_COMPILER_FLAG("-mtune=core2" XFLAGS_OPT_X86)
        if (XFLAGS_OPT_X86)
            set(GMXC_C_FLAGS "-mtune=core2 ${GMXC_C_FLAGS}")
        endif()
        CHECK_C_COMPILER_FLAG("-mtune=itanium2" XFLAGS_OPT_IA64)
        if (XFLAGS_OPT_IA64)
            set(GMXC_C_FLAGS "-mtune=itanium2 ${GMXC_C_FLAGS}")
        endif()
    else()
        CHECK_C_COMPILER_FLAG("/Qip " XFLAGS_OPT)
        if (XFLAGS_OPT)
            set(GMXC_C_FLAGS "/Qip ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_C_COMPILER_FLAG("/arch:SSE2" XFLAGS_OPT_SSE2)
        if (XFLAGS_OPT_SSE2)
            set(GMXC_C_FLAGS "/arch:SSE2 ${GMXC_C_FLAGS}")
        endif()
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if (NOT WIN32) 
        CHECK_CXX_COMPILER_FLAG("-ip -w -funroll-all-loops -std=gnu99" XXFLAGS_OPT)
        if (XXFLAGS_OPT)
            set(GMXC_CXX_FLAGS "-ip -w -funroll-all-loops ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_CXX_COMPILER_FLAG("-msse2" XXFLAGS_OPT_SSE2)
        if (XXFLAGS_OPT_SSE2)
            set(GMXC_CXX_FLAGS "-msse2 ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_CXX_COMPILER_FLAG("-mtune=pentium4" XXFLAGS_OPT_X86)
        if (XXFLAGS_OPT_X86)
            set(GMXC_CXX_FLAGS "-mtune=pentium4 ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_CXX_COMPILER_FLAG("-mtune=core2" XXFLAGS_OPT_X86)
        if (XXFLAGS_OPT_X86)
            set(GMXC_CXX_FLAGS "-mtune=core2 ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_CXX_COMPILER_FLAG("-mtune=itanium2" XXFLAGS_OPT_IA64)
        if (XXFLAGS_OPT_IA64)
            set(GMXC_CXX_FLAGS "-mtune=itanium2 ${GMXC_CXX_FLAGS}")
        endif()
    else()
        CHECK_CXX_COMPILER_FLAG("/Qip " XXFLAGS_OPT)
        if (XXFLAGS_OPT)
            set(GMXC_CXX_FLAGS "/Qip ${GMXC_CXX_FLAGS}")
        endif()
        CHECK_CXX_COMPILER_FLAG("/arch:SSE2" XXFLAGS_OPT_SSE2)
        if (XXFLAGS_OPT_SSE2)
            set(GMXC_CXX_FLAGS "/arch:SSE2 ${GMXC_CXX_FLAGS}")
        endif()
    endif()
endif()

#TODO: other compilers, like xlc

# now actually set the flags:
if (NOT DEFINED GMXCFLAGS_SET)
    set(GMXCFLAGS_SET true CACHE INTERNAL "Whether to reset the C flags" FORCE)
    # C
    set(CMAKE_C_FLAGS "${GMXC_C_FLAGS} ${CMAKE_C_FLAGS}" 
        CACHE STRING "Flags used by the compiler during all build types" FORCE)
    set(CMAKE_C_FLAGS_RELEASE "${GMXC_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS_RELEASE}" 
        CACHE STRING "Flags used by the compiler during release build" FORCE)
    # C++
    set(CMAKE_CXX_FLAGS "${GMXC_CXX_FLAGS} ${CMAKE_CXX_FLAGS}" 
        CACHE STRING "Flags used by the compiler during all build types" FORCE)
    set(CMAKE_CXX_FLAGS_RELEASE "${GMXC_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE}" 
        CACHE STRING "Flags used by the compiler during all release build" FORCE)
endif()

ENDMACRO(gmx_c_flags)
