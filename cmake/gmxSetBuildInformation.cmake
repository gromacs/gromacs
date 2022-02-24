#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2012- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.

# Check the username performing the build, as well as date, time, and
# build CPU features.
#
# The following variables will be set to the user/host/cpu used for
# configuration, or anonymous/unknown if it cannot be detected
# (Windows).
#
# BUILD_CPU_VENDOR
# BUILD_CPU_BRAND
# BUILD_CPU_FAMILY
# BUILD_CPU_MODEL
# BUILD_CPU_STEPPING
# BUILD_CPU_FEATURES
#

include(gmxDetectCpu)

function(gmx_set_build_information)
    # Set up some defaults that will usually be overwritten
    if(CMAKE_CROSSCOMPILING)
        set(_reason ", cross-compiled")
    endif()

    # Run the cpu detection. If it produces an empty output, set a
    # local value in the parent scope with a suitable fallback (which
    # hides the cached value).

    macro(gmx_get_build_cpu_string TYPE DEFAULT_VALUE)
        string(TOUPPER ${TYPE} UPPERTYPE)
        gmx_run_cpu_detection(${TYPE})
        set(OUTPUT_VALUE "${DEFAULT_VALUE}")
        if (CPU_DETECTION_${UPPERTYPE})
            set(OUTPUT_VALUE ${CPU_DETECTION_${UPPERTYPE}})
        endif()
        set(BUILD_CPU_${UPPERTYPE} ${OUTPUT_VALUE} PARENT_SCOPE)
    endmacro()

    gmx_get_build_cpu_string(vendor   "Unknown${_reason}")
    gmx_get_build_cpu_string(brand    "Unknown${_reason}")
    gmx_get_build_cpu_string(family   "0")
    gmx_get_build_cpu_string(model    "0")
    gmx_get_build_cpu_string(stepping "0")
    gmx_get_build_cpu_string(features "Unknown${_reason}")
endfunction()
