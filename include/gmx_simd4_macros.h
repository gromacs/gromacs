/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

/* The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code with a SIMD width of 4.
 * To support a new architecture, adding macros here should be all
 * that is needed.
 *
 * Note that this file is intended only for SIMD operations that require
 * a SIMD width of 4. In general gmx_simd_macros.h provides wider hardware
 * support, more functionality and higher performance, but the SIMD width is
 * not necessarily equal to 4.
 */

#ifdef _gmx_simd4_macros_h_
#error "gmx_simd4_macros.h included twice"
#else
#define _gmx_simd4_macros_h_


/* The SIMD width here is always 4, since that is the whole point */
#define GMX_SIMD4_WIDTH  4


#if defined GMX_SIMD4_SINGLE || defined GMX_SIMD4_DOUBLE
/* Precision set before inclusion, honour that request */
#else
/* Match precision to the Gromacs real precision */
#ifdef GMX_DOUBLE
#define GMX_SIMD4_DOUBLE
#else
#define GMX_SIMD4_SINGLE
#endif
#endif

#ifdef GMX_SIMD4_DOUBLE
typedef double  gmx_simd4_real;
#endif
#ifdef GMX_SIMD4_SINGLE
typedef float   gmx_simd4_real;
#endif

/* Uncomment the next line, without other SIMD active, for testing plain-C */
/* #define GMX_SIMD4_REFERENCE_PLAIN_C */
#ifdef GMX_SIMD4_REFERENCE_PLAIN_C
/* Plain C SIMD reference implementation, also serves as documentation */
#define GMX_HAVE_SIMD4_MACROS

/* Include plain-C reference implementation, also serves as documentation */
#include "gmx_simd4_ref.h"

/* float/double SIMD register type */
#define gmx_simd4_pr  gmx_simd4_ref_pr

/* boolean SIMD register type */
#define gmx_simd4_pb  gmx_simd4_ref_pb

#define gmx_simd4_load_pr       gmx_simd4_ref_load_pr
#define gmx_simd4_set1_pr       gmx_simd4_ref_set1_pr
#define gmx_simd4_setzero_pr    gmx_simd4_ref_setzero_pr
#define gmx_simd4_store_pr      gmx_simd4_ref_store_pr

/* Unaligned load+store are not required,
 * but they can speed up the PME spread+gather operations.
 */
#define GMX_SIMD4_HAVE_UNALIGNED
#ifdef GMX_SIMD4_HAVE_UNALIGNED
#define gmx_simd4_loadu_pr      gmx_simd4_ref_load_pr
#define gmx_simd4_storeu_pr     gmx_simd4_ref_store_pr
#endif

#define gmx_simd4_add_pr        gmx_simd4_ref_add_pr
#define gmx_simd4_sub_pr        gmx_simd4_ref_sub_pr
#define gmx_simd4_mul_pr        gmx_simd4_ref_mul_pr
/* For the FMA macros below, aim for c=d in code, so FMA3 uses 1 instruction */
#define gmx_simd4_madd_pr       gmx_simd4_ref_madd_pr
#define gmx_simd4_nmsub_pr      gmx_simd4_ref_nmsub_pr

#define gmx_simd4_dotproduct3   gmx_simd4_ref_dotproduct3

#define gmx_simd4_min_pr        gmx_simd4_ref_min_pr
#define gmx_simd4_max_pr        gmx_simd4_ref_max_pr

#define gmx_simd4_blendzero_pr  gmx_simd4_ref_blendzero_pr

/* Comparison */
#define gmx_simd4_cmplt_pr      gmx_simd4_ref_cmplt_pr

/* Logical operations on SIMD booleans */
#define gmx_simd4_and_pb        gmx_simd4_ref_and_pb
#define gmx_simd4_or_pb         gmx_simd4_ref_or_pb

/* Returns a single int (0/1) which tells if any of the 4 booleans is True */
#define gmx_simd4_anytrue_pb    gmx_simd4_ref_anytrue_pb

#endif /* GMX_SIMD4_REFERENCE_PLAIN_C */


/* The same SIMD macros can be translated to SIMD intrinsics (and compiled
 * to instructions for) different SIMD width and float precision.
 *
 * On x86: The gmx_simd4 prefix is replaced by _mm_ or _mm256_ (SSE or AVX).
 * The _pr suffix is replaced by _ps or _pd (for single or double precision).
 * Compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */


#ifdef GMX_X86_SSE2
/* This is for general x86 SIMD instruction sets that also support SSE2 */

#ifdef GMX_SIMD4_SINGLE
#define GMX_HAVE_SIMD4_MACROS
#endif

#ifdef GMX_SIMD4_DOUBLE
/* Note that here we will use 256-bit SIMD with GMX_X86_AVX_128_FMA.
 * This is inconsistent naming wise, but should give the best performance.
 */
#if defined GMX_X86_AVX_128_FMA || defined GMX_X86_AVX_256
#define GMX_HAVE_SIMD4_MACROS
#endif
#endif

#ifdef GMX_HAVE_SIMD4_MACROS

#if defined GMX_X86_AVX_128_FMA || defined GMX_X86_AVX_256

#include <immintrin.h>
#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h> /* FMA */
#endif
#ifdef HAVE_INTRIN_H
#include <intrin.h> /* FMA MSVC */
#endif

#else
#ifdef GMX_X86_SSE4_1
#include <smmintrin.h>
#else
/* We only have SSE2 */
#include <emmintrin.h>
#endif
#endif

#ifdef GMX_SIMD4_SINGLE

#define gmx_simd4_pr  __m128

#define gmx_simd4_pb  __m128

#define gmx_simd4_load_pr       _mm_load_ps
#define gmx_simd4_set1_pr       _mm_set1_ps
#define gmx_simd4_setzero_pr    _mm_setzero_ps
#define gmx_simd4_store_pr      _mm_store_ps

/* Some old AMD processors could have problems with unaligned loads+stores */
#ifndef GMX_FAHCORE
#define GMX_SIMD4_HAVE_UNALIGNED
#endif
#ifdef GMX_SIMD4_HAVE_UNALIGNED
#define gmx_simd4_loadu_pr      _mm_loadu_ps
#define gmx_simd4_storeu_pr     _mm_storeu_ps
#endif

#define gmx_simd4_add_pr        _mm_add_ps
#define gmx_simd4_sub_pr        _mm_sub_ps
#define gmx_simd4_mul_pr        _mm_mul_ps

#ifdef GMX_X86_AVX_128_FMA
#define gmx_simd4_madd_pr(a, b, c)   _mm_macc_ps(a, b, c)
#define gmx_simd4_nmsub_pr(a, b, c)  _mm_nmacc_ps(a, b, c)
#else
#define gmx_simd4_madd_pr(a, b, c)   _mm_add_ps(c, _mm_mul_ps(a, b))
#define gmx_simd4_nmsub_pr(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#endif

static inline float gmx_simd4_dotproduct3(__m128 a, __m128 b)
#ifdef GMX_X86_SSE4_1
{
    float dp;

    /* SSE4.1 dot product of components 0,1,2, stored in component 0 */
    _mm_store_ss(&dp, _mm_dp_ps(a, b, 0x71));

    return dp;
}
#else
{
    float        dp_array[7], *dp;

    /* Generate an aligned pointer */
    dp = (float *)(((size_t)(dp_array+3)) & (~((size_t)15)));

    _mm_store_ps(dp, _mm_mul_ps(a, b));

    return dp[0] + dp[1] + dp[2];
}
#endif

#define gmx_simd4_min_pr        _mm_min_ps
#define gmx_simd4_max_pr        _mm_max_ps

#define gmx_simd4_blendzero_pr  _mm_and_ps

#define gmx_simd4_cmplt_pr      _mm_cmplt_ps
#define gmx_simd4_and_pb        _mm_and_ps
#define gmx_simd4_or_pb         _mm_or_ps

#define gmx_simd4_anytrue_pb    _mm_movemask_ps

#endif /* GMX_SIMD4_SINGLE */


#ifdef GMX_SIMD4_DOUBLE

#define gmx_simd4_pr  __m256d

#define gmx_simd4_pb  __m256d

#define gmx_simd4_load_pr       _mm256_load_pd
#define gmx_simd4_set1_pr       _mm256_set1_pd
#define gmx_simd4_setzero_pr    _mm256_setzero_pd
#define gmx_simd4_store_pr      _mm256_store_pd

#define GMX_SIMD4_HAVE_UNALIGNED
#define gmx_simd4_loadu_pr      _mm256_loadu_pd
#define gmx_simd4_storeu_pr     _mm256_storeu_pd

#define gmx_simd4_add_pr        _mm256_add_pd
#define gmx_simd4_sub_pr        _mm256_sub_pd
#define gmx_simd4_mul_pr        _mm256_mul_pd
#ifdef GMX_X86_AVX_128_FMA
#define gmx_simd4_madd_pr(a, b, c)   _mm256_macc_pd(a, b, c)
#define gmx_simd4_nmsub_pr(a, b, c)  _mm256_nmacc_pd(a, b, c)
#else
#define gmx_simd4_madd_pr(a, b, c)   _mm256_add_pd(c, _mm256_mul_pd(a, b))
#define gmx_simd4_nmsub_pr(a, b, c)  _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#endif
#define gmx_simd4_min_pr        _mm256_min_pd
#define gmx_simd4_max_pr        _mm256_max_pd

#define gmx_simd4_blendzero_pr  _mm256_and_pd

/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_simd4_cmplt_pr(x, y) _mm256_cmp_pd(x, y, 0x11)
#define gmx_simd4_and_pb        _mm256_and_pd
#define gmx_simd4_or_pb         _mm256_or_pd

#define gmx_simd4_anytrue_pb    _mm256_movemask_pd

#endif /* GMX_SIMD4_DOUBLE */


#endif /* GMX_HAVE_SIMD4_MACROS */


#endif /* GMX_X86_SSE2 */


#ifdef GMX_HAVE_SIMD4_MACROS
/* Generic functions to extract a SIMD4 aligned pointer from a pointer x.
 * x should have at least GMX_SIMD4_WIDTH=4 elements extra compared
 * to how many you want to use, to avoid indexing outside the aligned region.
 */

static gmx_inline gmx_simd4_real *
gmx_simd4_align_real(const gmx_simd4_real *x)
{
    return (gmx_simd4_real *)(((size_t)((x)+GMX_SIMD4_WIDTH)) & (~((size_t)(GMX_SIMD4_WIDTH*sizeof(gmx_simd4_real)-1))));
}
#endif


#endif /* _gmx_simd4_macros_h_ */
