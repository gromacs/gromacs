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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_COMMON_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_COMMON_H

/* Please see documentation in gromacs/simd/simd.h for details. */

/* Inherit parts of AVX_128_FMA from SSE4.1 */
#include "gromacs/simd/impl_x86_sse4_1/impl_x86_sse4_1.h"

/* Override some capability definitions for things added in AVX over SSE4.1 */
#undef  GMX_SIMD_HAVE_FMA
#define GMX_SIMD_HAVE_FMA          1
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FRACTION     1
#undef  GMX_SIMD4_HAVE_DOUBLE
#define GMX_SIMD4_HAVE_DOUBLE      1 /* We can use 256-bit operations for this */

/* Work around gcc bug with wrong type for mask formal parameter to maskload/maskstore */
#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), _mm_castsi128_ps(mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)   _mm_maskstore_ps((mem), _mm_castsi128_ps(mask), (x))
#else
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), (mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)   _mm_maskstore_ps((mem), (mask), (x))
#endif

#endif  /* GMX_SIMD_IMPL_X86_AVX_128_FMA_COMMON_H */
