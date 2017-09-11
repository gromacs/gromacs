#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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

# Sources that should always be built
file(GLOB NONBONDED_SOURCES *.cpp nb_kernel_c/*.c)

if("${GMX_SIMD_ACTIVE}" STREQUAL "SSE2" AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_SSE2_SINGLE_SOURCES nb_kernel_sse2_single/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "SSE4.1" AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_SSE4_1_SINGLE_SOURCES nb_kernel_sse4_1_single/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_128_FMA" AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_AVX_128_FMA_SINGLE_SOURCES nb_kernel_avx_128_fma_single/*.c)
endif()

if((("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_256")
    OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX2_256")
    OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512")
    OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512_KNL"))
        AND NOT GMX_DOUBLE)
    file(GLOB NONBONDED_AVX_256_SINGLE_SOURCES nb_kernel_avx_256_single/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "SSE2" AND GMX_DOUBLE)
    file(GLOB NONBONDED_SSE2_DOUBLE_SOURCES nb_kernel_sse2_double/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "SSE4.1" AND GMX_DOUBLE)
    file(GLOB NONBONDED_SSE4_1_DOUBLE_SOURCES nb_kernel_sse4_1_double/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_128_FMA" AND GMX_DOUBLE)
    file(GLOB NONBONDED_AVX_128_FMA_DOUBLE_SOURCES nb_kernel_avx_128_fma_double/*.c)
endif()

if((("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_256")
            OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX2_256")
            OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512")
            OR ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512_KNL"))
        AND GMX_DOUBLE)
    file(GLOB NONBONDED_AVX_256_DOUBLE_SOURCES nb_kernel_avx_256_double/*.c)
endif()

if("${GMX_SIMD_ACTIVE}" STREQUAL "SPARC64_HPC_ACE" AND GMX_DOUBLE)
    file(GLOB NONBONDED_SPARC64_HPC_ACE_DOUBLE_SOURCES nb_kernel_sparc64_hpc_ace_double/*.c)
endif()


# These sources will be used in the parent directory's CMakeLists.txt
set(NONBONDED_SOURCES ${NONBONDED_SOURCES} ${NONBONDED_SSE2_SINGLE_SOURCES} ${NONBONDED_SSE4_1_SINGLE_SOURCES} ${NONBONDED_AVX_128_FMA_SINGLE_SOURCES} ${NONBONDED_AVX_256_SINGLE_SOURCES} ${NONBONDED_SSE2_DOUBLE_SOURCES} ${NONBONDED_SSE4_1_DOUBLE_SOURCES} ${NONBONDED_AVX_128_FMA_DOUBLE_SOURCES} ${NONBONDED_AVX_256_DOUBLE_SOURCES} ${NONBONDED_SPARC64_HPC_ACE_DOUBLE_SOURCES} PARENT_SCOPE)



