/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_COMMON_H
#define GMX_SIMD_IMPL_X86_AVX_512F_COMMON_H

#include <math.h>

#include <immintrin.h>

/* Intel AVX-512F */

/* This implementation works on the Intel compiler (version 14.0.2), but gcc-4.9
 * is still missing the cast operations needed to reinterpret between
 * 512-bit SIMD variables.
 * In addition, the truncate intrinsic is not yet available on gcc, and if
 * we try to implement it through _mm512_roundscale_round_ps() we get strange
 * errors about invalid immediate values when compiling Gromacs, but not a
 * simple test program. Hopefully these glitches will have been fixed in gcc
 * by the time AVX-512F-capable hardware appears.
 *
 * A couple of additional notes:
 *
 * - We avoid using the reduce calls provided by the intel compiler, since they
 *   are not instructions but software routines, and not implemented on gcc.
 * - gcc-4.9 does not yet provide int2mask and mask2int. This will be trivial
 *   to accomplish with plain casts, though.
 * - gcc-4.9 can do conversions between float/integer/double with simple casts
 *   like (__m512) as long as the variables are of the same size. However, it
 *   does not work to cast e.g. __m512 to __m128 this way. Since the sharing
 *   of xmm/ymm/zmm registers is more or less the whole idea with AVX-512
 *   (compared to MIC), this will have to be fixed in gcc.
 */

/* Capability definitions for AVX-512 SIMD. */
#define GMX_SIMD                             1
#define GMX_SIMD_HAVE_FLOAT                  1
#define GMX_SIMD_HAVE_DOUBLE                 1
#define GMX_SIMD_HAVE_LOADU                  1
#define GMX_SIMD_HAVE_STOREU                 1
#define GMX_SIMD_HAVE_LOGICAL                1
#define GMX_SIMD_HAVE_FMA                    1
#define GMX_SIMD_HAVE_FRACTION               0
#define GMX_SIMD_HAVE_FINT32                 1
/* Technically it is straightforward to emulate extract on AVX-512F through
 * memory operations, but when applied to 16 elements as part of a table lookup
 * it will be faster to just store the entire vector once, so we avoid setting it.
 */
#define GMX_SIMD_HAVE_FINT32_EXTRACT         0
#define GMX_SIMD_HAVE_FINT32_LOGICAL         1
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS     1
#define GMX_SIMD_HAVE_DINT32                 1
#define GMX_SIMD_HAVE_DINT32_EXTRACT         0
#define GMX_SIMD_HAVE_DINT32_LOGICAL         1
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS     1
#define GMX_SIMD4_HAVE_FLOAT                 1
#define GMX_SIMD4_HAVE_DOUBLE                1

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH                16
#define GMX_SIMD_DOUBLE_WIDTH                8
#define GMX_SIMD_FINT32_WIDTH               16
#define GMX_SIMD_DINT32_WIDTH                8
#define GMX_SIMD4_WIDTH                      4
#define GMX_SIMD_RSQRT_BITS                 14
#define GMX_SIMD_RCP_BITS                   14

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_COMMON_H */
