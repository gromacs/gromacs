#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013, by the GROMACS development team, led by
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

# Test C flags FLAGS, and set VARIABLE to true if they work. Also add the
# flags to C_FLAGS_VAR.
MACRO(GMX_TEST_CFLAG VARIABLE FLAGS C_FLAGS_VAR)
    IF(NOT DEFINED ${VARIABLE})
        CHECK_C_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF(NOT DEFINED ${VARIABLE})
    IF (${VARIABLE})
        SET (${C_FLAGS_VAR} "${FLAGS} ${${C_FLAGS_VAR}}")
    ENDIF (${VARIABLE})
ENDMACRO(GMX_TEST_CFLAG VARIABLE FLAGS C_FLAGS_VAR)

# Test C++ flags FLAGS, and set VARIABLE to true if they work. Also add the
# flags to CXX_FLAGS_VAR.
MACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXX_FLAGS_VAR)
    IF(NOT DEFINED ${VARIABLE} AND CMAKE_CXX_COMPILER_LOADED)
        CHECK_CXX_COMPILER_FLAG("${FLAGS}" ${VARIABLE})
    ENDIF(NOT DEFINED ${VARIABLE} AND CMAKE_CXX_COMPILER_LOADED)
    IF (${VARIABLE})
        SET (${CXX_FLAGS_VAR} "${FLAGS} ${${CXX_FLAGS_VAR}}")
    ENDIF (${VARIABLE})
ENDMACRO(GMX_TEST_CXXFLAG VARIABLE FLAGS CXX_FLAGS_VAR)

# TODO: the variable in the above macros that reports whether the
# flags that are generated actually work is never checked in the
# code below.

# This is the actual exported function to be called 
MACRO(gmx_manage_compiler_flags)

    include(CheckCCompilerFlag)
    include(CheckCXXCompilerFlag)

    # gcc
    if(CMAKE_COMPILER_IS_GNUCC)
        #flags are added in reverse order and -Wno* need to appear after -Wall
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(C_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_C_FLAGS)
        endif()
        GMX_TEST_CFLAG(C_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_C_FLAGS)
        GMX_TEST_CFLAG(C_FLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMX_DEFAULT_C_FLAGS)
        # new in gcc 4.5
        GMX_TEST_CFLAG(C_FLAGS_EXCESS_PREC "-fexcess-precision=fast" GMX_DEFAULT_C_FLAGS_RELEASE)
        GMX_TEST_CFLAG(C_FLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                       GMX_DEFAULT_C_FLAGS_RELEASE)
        GMX_TEST_CFLAG(C_FLAGS_NOINLINE "-fno-inline" GMX_DEFAULT_C_FLAGS_DEBUG)
    endif()
    # g++
    if(CMAKE_COMPILER_IS_GNUCXX)
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXX_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_CXX_FLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_CXX_FLAGS)
        GMX_TEST_CXXFLAG(CXX_FLAGS_WARN_EXTRA "-Wextra -Wno-missing-field-initializers -Wno-sign-compare" GMX_DEFAULT_CXX_FLAGS)
      # new in gcc 4.5
        GMX_TEST_CXXFLAG(CXX_FLAGS_EXCESS_PREC "-fexcess-precision=fast" GMX_DEFAULT_CXX_FLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXX_FLAGS_COPT "-fomit-frame-pointer -funroll-all-loops"
                         GMX_DEFAULT_CXX_FLAGS_RELEASE)
        GMX_TEST_CXXFLAG(CXX_FLAGS_NOINLINE "-fno-inline" GMX_DEFAULT_CXX_FLAGS_DEBUG)
    endif()

    # icc
    if (CMAKE_C_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CFLAG(C_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_C_FLAGS)
            endif()
            GMX_TEST_CFLAG(C_FLAGS_WARN "-Wall" GMX_DEFAULT_C_FLAGS)
            GMX_TEST_CFLAG(C_FLAGS_STDGNU "-std=gnu99" GMX_DEFAULT_C_FLAGS)
            GMX_TEST_CFLAG(C_FLAGS_OPT "-ip -funroll-all-loops" GMX_DEFAULT_C_FLAGS_RELEASE)
        else()
            GMX_TEST_CFLAG(C_FLAGS_WARN "/W2" GMX_DEFAULT_C_FLAGS)
            GMX_TEST_CFLAG(C_FLAGS_X86 "/Qip" GMX_DEFAULT_C_FLAGS_RELEASE)
        endif()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        if (NOT WIN32) 
            if(NOT GMX_OPENMP)
                GMX_TEST_CXXFLAG(CXX_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_CXX_FLAGS)
            endif()
            GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "-Wall" GMX_DEFAULT_CXX_FLAGS)
            GMX_TEST_CXXFLAG(CXX_FLAGS_OPT "-ip -funroll-all-loops" GMX_DEFAULT_CXX_FLAGS_RELEASE)
        else()
            GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "/W2" GMX_DEFAULT_CXX_FLAGS)
            GMX_TEST_CXXFLAG(CXX_FLAGS_X86 "/Qip" GMX_DEFAULT_CXX_FLAGS_RELEASE)
        endif()
    endif()

    # pgi
    if (CMAKE_C_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CFLAG(C_FLAGS_OPT "-fastsse" GMX_DEFAULT_C_FLAGS_RELEASE)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        GMX_TEST_CXXFLAG(CXX_FLAGS_OPT "-fastsse" GMX_DEFAULT_CXX_FLAGS_RELEASE)
    endif()

    # Pathscale
    if (CMAKE_C_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(C_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_C_FLAGS)
        endif()
        GMX_TEST_CFLAG(C_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_C_FLAGS)
        GMX_TEST_CFLAG(C_FLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMX_DEFAULT_C_FLAGS_RELEASE)
        GMX_TEST_CFLAG(C_FLAGS_LANG "-std=gnu99" GMX_DEFAULT_C_FLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "PathScale")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXX_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_CXX_FLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_CXX_FLAGS)
        GMX_TEST_CXXFLAG(CXX_FLAGS_OPT "-OPT:Ofast -fno-math-errno -ffast-math" 
                         GMX_DEFAULT_CXX_FLAGS_RELEASE)
    endif()

    # xlc
    if (CMAKE_C_COMPILER_ID MATCHES "XL")
        GMX_TEST_CFLAG(C_FLAGS_OPT "-qarch=auto -qtune=auto" GMX_DEFAULT_C_FLAGS)
        GMX_TEST_CFLAG(C_FLAGS_LANG "-qlanglvl=extc99" GMX_DEFAULT_C_FLAGS)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        GMX_TEST_CXXFLAG(CXX_FLAGS_OPT "-qarch=auto -qtune=auto" GMX_DEFAULT_CXX_FLAGS)
    endif()

    #msvc
    if (MSVC)
        # disable warnings for: 
        #      inconsistent dll linkage
        #      forcing value to bool (for C++)
        GMX_TEST_CFLAG(C_FLAGS_WARN "/wd4273" GMX_DEFAULT_C_FLAGS)
        GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "/wd4273 /wd4800" GMX_DEFAULT_CXX_FLAGS)
    endif()

    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CFLAG(C_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_C_FLAGS)
        endif()
        GMX_TEST_CFLAG(C_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_C_FLAGS)
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(NOT GMX_OPENMP)
            GMX_TEST_CXXFLAG(CXX_FLAGS_PRAGMA "-Wno-unknown-pragmas" GMX_DEFAULT_CXX_FLAGS)
        endif()
        GMX_TEST_CXXFLAG(CXX_FLAGS_WARN "-Wall -Wno-unused -Wunused-value" GMX_DEFAULT_CXX_FLAGS)
    endif()

    # now actually set the flags:
    # C
    if ( NOT GMX_SKIP_DEFAULT_FLAGS )
        set(CMAKE_C_FLAGS
            "${GMX_DEFAULT_C_FLAGS} ${USER_APPEND_C_FLAGS} ${CMAKE_C_FLAGS}")
        set(CMAKE_C_FLAGS_RELEASE
            "${GMX_DEFAULT_C_FLAGS_RELEASE} ${USER_APPEND_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS_RELEASE}")
        set(CMAKE_C_FLAGS_DEBUG
            "${GMX_DEFAULT_C_FLAGS_DEBUG} ${USER_APPEND_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_DEBUG}")
    endif()

    # C++
    if ( NOT GMX_SKIP_DEFAULT_FLAGS)
        set(CMAKE_CXX_FLAGS
            "${GMX_DEFAULT_CXX_FLAGS} ${USER_APPEND_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS_RELEASE 
            "${GMX_DEFAULT_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE}")
        set(CMAKE_CXX_FLAGS_DEBUG 
            "${GMX_DEFAULT_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_DEBUG}")
    endif()

    # Add advanced cache entries for the user flags, whether the user
    # has set them or not, so that they're accessible in ccmake.
    set(USER_APPEND_C_FLAGS ${USER_APPEND_C_FLAGS} CACHE STRING "Extra C compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for all build types.")
    mark_as_advanced(USER_APPEND_C_FLAGS_RELEASE)
    set(USER_APPEND_C_FLAGS_RELEASE ${USER_APPEND_C_FLAGS_RELEASE} CACHE STRING "Extra C compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for Release build types.")
    mark_as_advanced(USER_APPEND_C_FLAGS_RELEASE)
    set(USER_APPEND_C_FLAGS_DEBUG ${USER_APPEND_C_FLAGS_DEBUG} CACHE STRING "Extra C compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for Debug build types.")
    mark_as_advanced(USER_APPEND_C_FLAGS_DEBUG)

    set(USER_APPEND_CXX_FLAGS ${USER_APPEND_CXX_FLAGS} CACHE STRING "Extra C++ compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for all build types.")
    mark_as_advanced(USER_APPEND_CXX_FLAGS_RELEASE)
    set(USER_APPEND_CXX_FLAGS_RELEASE ${USER_APPEND_CXX_FLAGS_RELEASE} CACHE STRING "Extra C++ compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for Release build types.")
    mark_as_advanced(USER_APPEND_CXX_FLAGS_RELEASE)
    set(USER_APPEND_CXX_FLAGS_DEBUG ${USER_APPEND_CXX_FLAGS_DEBUG} CACHE STRING "Extra C++ compiler flags provided by you that will be appended to the flags generated by GROMACS, and used for Debug build types.")
    mark_as_advanced(USER_APPEND_CXX_FLAGS_DEBUG)

ENDMACRO()

