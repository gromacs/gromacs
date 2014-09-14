#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

# Due to a bug, gcc 4.4.x crashes when compiling listed-forces/bonded.cpp with -O3 and
# -fopenmp, but strangely it does not crash with -O2 + all additional options.
# -O3 uses. Therefore, for the affected files, when compiling in release mode,
# we override -O3 with -O2 and add the additional option.
#

# Considering compiler version and build configuration, check if the workaround
# is needed to avoid gcc crash.
macro(gmx_check_gcc44_bug_workaround_needed OUT_VAR)
    if(CMAKE_COMPILER_IS_GNUCC AND
    CMAKE_C_COMPILER_VERSION VERSION_GREATER "4.3.999" AND CMAKE_C_COMPILER_VERSION VERSION_LESS "4.4.999")

        set(_gcc44_workaround FALSE)

        # only apply the workaround if we are actually using -O3
        string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
        if ("${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${_build_type}}" MATCHES ".*-O3.*" AND
            GMX_OPENMP)
            if(GMX_DISABLE_GCC44_BUG_WORKAROUND)
                set(_msg "gcc ${CMAKE_C_COMPILER_VERSION} detected, using -O3, but workaround for optimization bug is disabled")
            else()
                set(_msg "gcc ${CMAKE_C_COMPILER_VERSION} detected, using -O3, will apply workaround for optimization bug (disable with GMX_DISABLE_GCC44_BUG_WORKAROUND)")
                set(_gcc44_workaround TRUE)
            endif()
            # only issues message if the value has changed
            if((NOT _gcc44_workaround AND ${OUT_VAR}) OR (_gcc44_workaround AND NOT ${OUT_VAR}))
                message(STATUS "${_msg}")
            endif()
        endif()

        set(${OUT_VAR} ${_gcc44_workaround} CACHE INTERNAL "Use gcc 4.4.x O3 optimization bug workaround" FORCE)
    endif()
endmacro()

# Apply workaround on the specified source file.
#
# This workaround does not seem to affect the performance in a measurable way.
macro(gmx_apply_gcc44_bug_workaround FILE_NAME)
    set_source_files_properties(
        ${FILE_NAME}
        PROPERTIES
        COMPILE_FLAGS "-O2 -finline-functions -funswitch-loops -fpredictive-commoning -fgcse-after-reload -ftree-vectorize -fipa-cp-clone"
        )
endmacro()
