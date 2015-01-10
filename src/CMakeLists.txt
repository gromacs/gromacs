#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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

######################################
# Output compiler and CFLAGS used
######################################
include(GetCompilerInfo.cmake)
get_compiler_info(C BUILD_C_COMPILER BUILD_CFLAGS)
get_compiler_info(CXX BUILD_CXX_COMPILER BUILD_CXXFLAGS)
if(GMX_GPU AND NOT GMX_USE_OPENCL)
    get_cuda_compiler_info(CUDA_NVCC_COMPILER_INFO CUDA_NVCC_COMPILER_FLAGS)
endif()

configure_file(config.h.cmakein config.h)
configure_file(gmxpre-config.h.cmakein gmxpre-config.h)
configure_file(buildinfo.h.cmakein buildinfo.h ESCAPE_QUOTES)

if (BUILD_TESTING)
    if(NOT GMX_DEVELOPER_BUILD)
        set(UNITTEST_TARGET_OPTIONS EXCLUDE_FROM_ALL)
    endif()
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(external/gmock-1.7.0)
    endif()
    include(testutils/TestMacros.cmake)
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(testutils)
    else()
        add_custom_target(unittests-notice
            ${CMAKE_COMMAND} -E echo "NOTE: Unit tests have not been run. You need to set GMX_BUILD_UNITTESTS=ON if you want to build and run them."
            DEPENDS run-ctest
            COMMENT "Unit tests disabled" VERBATIM)
        add_dependencies(check unittests-notice)
    endif()
endif()

add_subdirectory(gromacs)
add_subdirectory(programs)

if (NOT GMX_FAHCORE)
    add_subdirectory(contrib)
endif()
