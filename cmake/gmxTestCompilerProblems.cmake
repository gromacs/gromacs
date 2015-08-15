#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
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

# Macro that runs through a number of tests for buggy compiler
# versions, or other potential problems.
macro(gmx_test_compiler_problems)

    # Warn if C and C++ compilers do not match
    if(NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_CXX_COMPILER_ID)
        message(WARNING "The ids of the C and C++ compilers do not match (${CMAKE_C_COMPILER_ID} and ${CMAKE_CXX_COMPILER_ID}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()
    if(NOT CMAKE_C_COMPILER_VERSION STREQUAL CMAKE_CXX_COMPILER_VERSION)
        message(WARNING "The versions of the C and C++ compilers do not match (${CMAKE_C_COMPILER_VERSION} and ${CMAKE_CXX_COMPILER_VERSION}, respectively). Mixing different C/C++ compilers can cause problems.")
    endif()

    # gcc 4.4.x is buggy and crashes when compiling some files with O3 and OpenMP on.
    # Detect here whether applying a workaround is needed and will apply it later
    # on the affected files. This test must come after gmx_c_flags(), since we
    # only want to enable the workaround when using the -O3 flag.
    include(gmxGCC44O3BugWorkaround)
    gmx_check_gcc44_bug_workaround_needed(GMX_USE_GCC44_BUG_WORKAROUND)

    # clang 3.0 is buggy for some unknown reason detected during adding
    # the SSE2 group kernels for GROMACS 4.6. If we ever work out what
    # that is, we should replace these tests with a compiler feature test,
    # update GROMACS Redmine task #1039 and perhaps report a clang bug.
    #
    # In the meantime, until we require CMake 2.8.10 we cannot rely on it to detect
    # the compiler version for us. So we need a manual check for clang 3.0.
    include(gmxDetectClang30)
    gmx_detect_clang_3_0(COMPILER_IS_CLANG_3_0)
    if(COMPILER_IS_CLANG_3_0)
        message(FATAL_ERROR "Your compiler is clang version 3.0, which is known to be buggy for GROMACS. Use a different compiler.")
    endif()

    # clang <=3.2 contains a bug that causes incorrect code to be generated for the
    # vfmaddps instruction and therefore the bug is triggered with AVX_128_FMA.
    # (see: http://llvm.org/bugs/show_bug.cgi?id=15040).
    # We can work around this by not using the integrated assembler (except on OS X
    # which has an outdated assembler that does not support AVX instructions).
    if (CMAKE_C_COMPILER_ID MATCHES "Clang" AND CMAKE_C_COMPILER_VERSION VERSION_LESS "3.3")
        set(GMX_USE_CLANG_C_FMA_BUG_WORKAROUND TRUE)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.3")
        set(GMX_USE_CLANG_CXX_FMA_BUG_WORKAROUND TRUE)
    endif()

    if (CMAKE_C_COMPILER_ID STREQUAL "PGI")
        message(WARNING "Currently tested PGI compiler versions (up to 15.7) generate binaries that do not pass all regression test, and the generated binaries are significantly slower than with GCC, ICC or Clang. For now we do not recommend PGI beyond development testing - make sure to run the regressiontests.")
    endif()

    if(CMAKE_COMPILER_IS_GNUCC AND
            (CMAKE_C_COMPILER_VERSION VERSION_LESS "4.9.0" OR CMAKE_SIZEOF_VOID_P EQUAL 8)
            AND (WIN32 OR CYGWIN)
            AND GMX_SIMD MATCHES "AVX" AND NOT GMX_SIMD STREQUAL AVX_128_FMA)
        message(WARNING "GCC on Windows (GCC older than 4.9 or any version when compiling for 64bit) with AVX (other than AVX_128_FMA) crashes. Choose a different GMX_SIMD or a different compiler.") # GCC bug 49001, 54412.
    endif()

    if(CMAKE_C_COMPILER_ID MATCHES "Clang" AND WIN32)
        if(CMAKE_VERSION VERSION_LESS 3.0.0)
            message(WARNING "Clang on Windows requires cmake 3.0.0")
        endif()
        if(CMAKE_C_COMPILER_VERSION VERSION_LESS 3.5.0)
            message(WARNING "Clang on Windows requires clang 3.5.0")
        endif()
    endif()

    if(CMAKE_C_COMPILER_ID MATCHES "Intel" AND CMAKE_C_COMPILER_VERSION VERSION_LESS "12.0.0")
        message(WARNING "Intel compilers before 12.0.0 are not routinely tested, so there may be problems. Version 11.1 with SSE4.1 is known to produce incorrect results. It is highly recommended to use a more up-to-date compiler. If you choose to use this version, make sure you run the regressiontests.")
    endif()

endmacro(gmx_test_compiler_problems)
