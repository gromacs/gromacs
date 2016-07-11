#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

# - Decide whether to use CPU cycle counters
#
# gmx_manage_cycle_counters()
#
# By default, we enable GMX_CYCLECOUNTERS for all architectures except ARMv7.
# On ARMv7, we enable it if we are not cross-compiling and can run a small
# test to confirm that the support is present in the kernel, otherwise we
# disable it.
#
macro(gmx_manage_cycle_counters)

    if(NOT DEFINED GMX_CYCLECOUNTERS)

        if(GMX_TARGET_ARMV7)

            if(NOT CMAKE_CROSSCOMPILING)

                try_run(ARMV7_COUNTER_RUN_VAR ARMV7_COUNTER_COMPILE_VAR
                        ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake/TestARMv7CycleCounters.cpp")

                # Enable cycle counter usage if the test ran fine and exited with 0 return code
                if(${ARMV7_COUNTER_COMPILE_VAR} AND ("${ARMV7_COUNTER_RUN_VAR}" EQUAL "0"))
                    set(GMX_CYCLECOUNTERS ON CACHE BOOL "Use CPU cycle counters timing")
                else()
                    set(GMX_CYCLECOUNTERS OFF CACHE BOOL "Use CPU cycle counters for timing")
                endif()

            else()

                # Disable cycle counters when cross-compiling for ARMv7
                set(GMX_CYCLECOUNTERS OFF CACHE BOOL "Use CPU cycle counters for timing")

            endif()

        else()

            # For now we (try to) enable cycle counters on all other platforms
            set(GMX_CYCLECOUNTERS ON CACHE BOOL "Use CPU cycle counters timing")

        endif()

        mark_as_advanced(GMX_CYCLECOUNTERS)

    endif()

endmacro()

