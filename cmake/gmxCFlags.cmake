MACRO(gmx_c_flags)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

if(CMAKE_COMPILER_IS_GNUCC)
# this one is set by cmake depending on the build targeb
#  CHECK_C_COMPILER_FLAG( "-O3" XFLAGS_O3)
#  if (XFLAGS_O3)
#    set(GMXC_C_FLAGS " -O3 ${GMXC_C_FLAGS}")
#  endif(XFLAGS_O3)

  CHECK_C_COMPILER_FLAG( "-Wall -Wno-unused" XFLAGS_WARN)
  if (XFLAGS_WARN)
    set(GMXC_C_FLAGS " -Wall -Wno-unused ${GMXC_C_FLAGS}")
  endif(XFLAGS_WARN)

  CHECK_C_COMPILER_FLAG( "-std=gnu99" XFLAGS_GNU99)
  if (XFLAGS_GNU99)
    set(GMXC_C_FLAGS "-std=gnu99 ${GMXC_C_FLAGS}")
  endif(XFLAGS_GNU99)

  CHECK_C_COMPILER_FLAG( "-march=native" XFLAGS_MARCH)
  if (XFLAGS_MARCH)
    set(GMXC_C_FLAGS "-march=native ${GMXC_C_FLAGS}")
  endif(XFLAGS_MARCH)
  # new in gcc 4.5
  CHECK_C_COMPILER_FLAG( "-fexcess-precision=fast" XFLAGS_EXCESS_PRECISION)
  if (XFLAGS_EXCESS_PRECISION)
    set(GMXC_C_FLAGS "-fexcess-precision=fast ${GMXC_C_FLAGS}")
  endif (XFLAGS_EXCESS_PRECISION)

  CHECK_C_COMPILER_FLAG( "-fomit-frame-pointer -finline-functions -funroll-all-loops" XFLAGS_OPT)
  if (XFLAGS_OPT)
    set(GMXC_C_FLAGS_RELEASE "-fomit-frame-pointer -finline-functions -funroll-all-loops ${GMXC_C_FLAGS_RELEASE}")
  endif(XFLAGS_OPT)
endif( CMAKE_COMPILER_IS_GNUCC )


if(CMAKE_COMPILER_IS_GNUCXX )
# this one is set by cmake depending on the build targeb
#  CHECK_CXX_COMPILER_FLAG( "-O3" XXFLAGS_O3)
#  if (XXFLAGS_O3)
##    set(GMXC_CXX_FLAGS " -O3 ${GMXC_CXX_FLAGS}")
#  endif(XXFLAGS_O3)

  CHECK_CXX_COMPILER_FLAG( "-Wall -Wno-unused" XXFLAGS_WARN)
  if (XXFLAGS_WARN)
    set(GMXC_CXX_FLAGS " -Wall -Wno-unused ${GMXC_CXX_FLAGS}")
  endif(XXFLAGS_WARN)

  CHECK_CXX_COMPILER_FLAG( "-std=gnu99" XXFLAGS_GNU99)
  if (XXFLAGS_GNU99)
    set(GMXC_CXX_FLAGS "-std=gnu99 ${GMXC_CXX_FLAGS}")
  endif(XXFLAGS_GNU99)

  CHECK_CXX_COMPILER_FLAG( "-march=native" XXFLAGS_MARCH)
  if (XXFLAGS_MARCH)
    set(GMXC_CXX_FLAGS "-march=native ${GMXC_CXX_FLAGS}")
  endif(XXFLAGS_MARCH)
  # new in gcc 4.5
  CHECK_CXX_COMPILER_FLAG( "-fexcess-precision=fast" XXFLAGS_EXCESS_PRECISION)
  if (XXFLAGS_EXCESS_PRECISION)
    set(GMXC_CXX_FLAGS "-fexcess-precision=fast ${GMXC_CXX_FLAGS}")
  endif (XXFLAGS_EXCESS_PRECISION)

  CHECK_CXX_COMPILER_FLAG( "-fomit-frame-pointer -finline-functions -funroll-all-loops" XXFLAGS_OPT)
  if (XXFLAGS_OPT)
    set(GMXC_CXX_FLAGS_RELEASE "-fomit-frame-pointer -finline-functions -funroll-all-loops ${GMXC_CXX_FLAGS_RELEASE}")
  endif(XXFLAGS_OPT)
endif(CMAKE_COMPILER_IS_GNUCXX )



if (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if (NOT WIN32) 
        CHECK_C_COMPILER_FLAG("-tpp7 -axW -ip -w -msse2 -funroll-all-loops -std=gnu99" XFLAGS_OPT)
        if (XFLAGS_OPT)
            set(GMXC_C_FLAGS "-tpp7 -axW -ip -w -msse2 -funroll-all-loops -std=gnu99 ${GMXC_C_FLAGS}")
        endif(XFLAGS_OPT)
    endif ()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    if (NOT WIN32) 
        CHECK_CXX_COMPILER_FLAG("-tpp7 -axW -ip -w -msse2 -funroll-all-loops -std=gnu99" XXFLAGS_OPT)
        if (XXFLAGS_OPT)
            set(GMXC_CXX_FLAGS "-tpp7 -axW -ip -w -msse2 -funroll-all-loops -std=gnu99 ${GMXC_CXX_FLAGS}")
        endif(XXFLAGS_OPT)
    endif ()
endif()


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
endif (NOT DEFINED GMXCFLAGS_SET)



ENDMACRO(gmx_c_flags)
