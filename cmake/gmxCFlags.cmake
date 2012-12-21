#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.
#

# Test C flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CFLAGSVAR.
MACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_C_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF(NOT DEFINED ${VARIABLE})
    IF (${VARIABLE})
        SET (${CFLAGSVAR} "${FLAGS} ${${CFLAGSVAR}}")
    ENDIF (${VARIABLE})
ENDMACRO(GMX_TEST_CFLAG VARIABLE FLAGS CFLAGSVAR)

# Test C++ flags FLAGS, and set VARIABLE to true if the work. Also add the
# flags to CXXFLAGSVAR.
MACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)
    IF(NOT DEFINED ${VARIABLE} AND CMAKE_CXX_COMPILER_LOADED)
        CHECK_CXX_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF(NOT DEFINED ${VARIABLE} AND CMAKE_CXX_COMPILER_LOADED)
    IF (${VARIABLE})
        SET (${CXXFLAGSVAR} "${FLAGS} ${${CXXFLAGSVAR}}")
    ENDIF (${VARIABLE})
ENDMACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXXFLAGSVAR)


# This is the actual exported function to be called 
MACRO(gmx_c_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # gcc
    if(CMAKE_COMPILER_IS_GNUCC)

        #Fix for LLVM OpenMP bug (redmine 900). Needs to run before OpenMP flags are set below.
        if(GMX_OPENMP)
            exec_program(${CMAKE_C_COMPILER} ARGS --version OUTPUT_VARIABLE _compiler_output)
            if(_compiler_output MATCHES "llvm.*4\\.2")
                message(STATUS "OpenMP multithreading not supported with llvm-gcc 4.2, disabled")
                set(GMX_OPENMP OFF CACHE BOOL
                    "OpenMP multithreading not not supported with llvm-gcc 4.2, disabled!" FORCE)
            endif()
        endif()

        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMXC_CFLAGS)
        # new in gcc 4.5
        GMX_TEST_CFLAG(CFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                       GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_NOINLINE "-fno-inline" GMXC_CFLAGS_DEBUG)
    endif()
    # g++
    if(CMAKE_COMPILER_IS_GNUCXX)
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMXC_CXXFLAGS)
      # new in gcc 4.5
        GMX_TEST_CXXFLAG(CXXFLAGS_EXCESS_PREC "-fexcess-precision=fast" GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                         GMXC_CXXFLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXXFLAGS_NOINLINE "-fno-inline" GMXC_CXXFLAGS_DEBUG)
    endif()

    # icc
    if (CMAKE_C_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
            endif()
            GMX_TEST_CFLAG(CFLAGS_WARN "-Wall" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_STDGNU "-std=gnu99" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_OPT "-ip -funroll-all-loops" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_SSE2 "-msse2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_X86 "-mtune=core2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_IA64 "-mtune=itanium2" GMXC_CFLAGS_RELEASE)
        else()
            #ignore INCREMENTAL:YES dll linkage warnings
            GMX_TEST_CFLAG(CFLAGS_WARN "/W2 /wd4224" GMXC_CFLAGS)
            GMX_TEST_CFLAG(CFLAGS_SSE2 "/arch:SSE2" GMXC_CFLAGS_RELEASE)
            GMX_TEST_CFLAG(CFLAGS_X86 "/Qip" GMXC_CFLAGS_RELEASE)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-ip -funroll-all-loops" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "-msse2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_X86 "-mtune=core2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_IA64 "-mtune=itanium2" 
                              GMXC_CXXFLAGS_RELEASE)
        else()
            GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/W2 /wd4224" GMXC_CXXFLAGS)
            GMX_TEST_CXXFLAG(CXXFLAGS_SSE2 "/arch:SSE2" GMXC_CXXFLAGS_RELEASE)
            GMX_TEST_CXXFLAG(CXXFLAGS_X86 "/Qip" GMXC_CXXFLAGS_RELEASE)
        endif()
    endif()

    # pgi
    if (CMAKE_C_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CFLAG(CFLAGS_OPT "-fastsse" GMXC_CFLAGS_RELEASE)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-fastsse" GMXC_CXXFLAGS_RELEASE)
    endif()

    # Pathscale
    if (CMAKE_C_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMXC_CFLAGS_RELEASE)
        GMX_TEST_CFLAG(CFLAGS_LANG "-std=gnu99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMXC_CXXFLAGS_RELEASE)
    endif()

    # xlc
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(CFLAGS_OPT "-qarch=auto -qtune=auto" GMXC_CFLAGS)
        GMX_TEST_CFLAG(CFLAGS_LANG "-qlanglvl=extc99" GMXC_CFLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        GMX_TEST_CXXFLAG(CXXFLAGS_OPT "-qarch=auto -qtune=auto" GMXC_CXXFLAGS)
    endif()

    #msvc
    if (MSVC)
        # disable warnings for: 
        #      inconsistent dll linkage
        GMX_TEST_CFLAG(CFLAGS_WARN "/wd4273" GMXC_CFLAGS)
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "/wd4273" GMXC_CXXFLAGS)
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(CFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CFLAGS)
        endif()
        GMX_TEST_CFLAG(CFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CFLAGS)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXXFLAGS_PRAGMA "-Wno-unknown-pragmas" GMXC_CXXFLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXXFLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMXC_CXXFLAGS)
    endif()

    # now actually set the flags:
    # C
    if ( NOT DEFINED GMXCFLAGS_SET AND NOT DEFINED ENV{CFLAGS} )
        set(GMXCFLAGS_SET true CACHE INTERNAL "Whether to reset the C flags" 
            FORCE)
        
        set(CMAKE_C_FLAGS "${GMXC_CFLAGS} ${CMAKE_C_FLAGS}" 
            CACHE STRING "Flags used by the compiler during all build types." 
            FORCE)
        set(CMAKE_C_FLAGS_RELEASE "${GMXC_CFLAGS_RELEASE} ${CMAKE_C_FLAGS_RELEASE}" 
            CACHE STRING "Flags used by the compiler during release builds." 
            FORCE)
        set(CMAKE_C_FLAGS_DEBUG "${GMXC_CFLAGS_DEBUG} ${CMAKE_C_FLAGS_DEBUG}" 
            CACHE STRING "Flags used by the compiler during debug builds." 
            FORCE)
    endif()

    # C++
    if ( NOT DEFINED GMXCXXFLAGS_SET AND NOT DEFINED ENV{CXXFLAGS} AND CMAKE_CXX_COMPILER_LOADED)
        set(GMXCXXFLAGS_SET true CACHE INTERNAL "Whether to reset the C++ flags" 
            FORCE)
        set(CMAKE_CXX_FLAGS "${GMXC_CXXFLAGS} ${CMAKE_CXX_FLAGS}" 
            CACHE STRING "Flags used by the compiler during all build types." 
            FORCE)
        set(CMAKE_CXX_FLAGS_RELEASE 
            "${GMXC_CXXFLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE}" 
            CACHE STRING "Flags used by the compiler during release builds." 
            FORCE)
        set(CMAKE_CXX_FLAGS_DEBUG 
            "${GMXC_CXXFLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG}" 
            CACHE STRING "Flags used by the compiler during debug builds." 
            FORCE)
    endif()
ENDMACRO(gmx_c_flags)

