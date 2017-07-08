/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_SIMD_IMPL_X86_AVX_256_DEFINITIONS_H
#define GMX_SIMD_IMPL_X86_AVX_256_DEFINITIONS_H

// Capability definitions for 256-bit AVX
#define GMX_SIMD                                1
#define GMX_SIMD_HAVE_FLOAT                     1
#define GMX_SIMD_HAVE_DOUBLE                    1
#define GMX_SIMD_HAVE_LOADU                     1
#define GMX_SIMD_HAVE_STOREU                    1
#define GMX_SIMD_HAVE_LOGICAL                   1
#define GMX_SIMD_HAVE_FMA                       0
#define GMX_SIMD_HAVE_FINT32_EXTRACT            1 // Emulated
#define GMX_SIMD_HAVE_FINT32_LOGICAL            0 // AVX1 cannot do 256-bit int shifts
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS        0 // AVX1 cannot do 256-bit int +,-,*
#define GMX_SIMD_HAVE_DINT32_EXTRACT            1 // Native, since we use __m128i
#define GMX_SIMD_HAVE_DINT32_LOGICAL            1
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS        1
#define GMX_SIMD_HAVE_NATIVE_COPYSIGN_FLOAT     0
#define GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_FLOAT   0
#define GMX_SIMD_HAVE_NATIVE_RCP_ITER_FLOAT     0
#define GMX_SIMD_HAVE_NATIVE_LOG_FLOAT          0
#define GMX_SIMD_HAVE_NATIVE_EXP2_FLOAT         0
#define GMX_SIMD_HAVE_NATIVE_EXP_FLOAT          0
#define GMX_SIMD_HAVE_NATIVE_COPYSIGN_DOUBLE    0
#define GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_DOUBLE  0
#define GMX_SIMD_HAVE_NATIVE_RCP_ITER_DOUBLE    0
#define GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE         0
#define GMX_SIMD_HAVE_NATIVE_EXP2_DOUBLE        0
#define GMX_SIMD_HAVE_NATIVE_EXP_DOUBLE         0
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT   1
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE  1
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT          1
#define GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE         0 // Not needed for width 4
#define GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT         1

#define GMX_SIMD4_HAVE_FLOAT                    1
#define GMX_SIMD4_HAVE_DOUBLE                   1

// Implementation details
#define GMX_SIMD_FLOAT_WIDTH                    8
#define GMX_SIMD_DOUBLE_WIDTH                   4
#define GMX_SIMD_FINT32_WIDTH                   8
#define GMX_SIMD_DINT32_WIDTH                   4
#define GMX_SIMD4_WIDTH                         4
#define GMX_SIMD_RSQRT_BITS                    11
#define GMX_SIMD_RCP_BITS                      11

#endif // GMX_SIMD_IMPL_X86_AVX_256_DEFINITIONS_H
