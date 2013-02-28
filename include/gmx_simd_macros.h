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
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding macros here should be (nearly)
 * all that is needed.
 */

/* Undefine all defines used below so we can include this file multiple times
 * with different settings from the same source file.
 */

/* NOTE: SSE2 acceleration does not include floor or blendv */

#undef GMX_SIMD_WIDTH_HERE

#undef gmx_epi32

/* float/double SIMD register type */
#undef gmx_mm_pr

#undef gmx_load_pr
#undef gmx_load1_pr
#undef gmx_set1_pr
#undef gmx_setzero_pr
#undef gmx_store_pr
/* Only used for debugging */
#undef gmx_storeu_pr

#undef gmx_add_pr
#undef gmx_sub_pr
#undef gmx_mul_pr
/* For the FMA macros below, aim for c=d in code, so FMA3 uses 1 instruction */
/* d = gmx_madd_pr(a,b,c): d = a*b + c, could use FMA3 or FMA4 */
#undef gmx_madd_pr
/* d = gmx_nmsub_pr(a,b,c): d = -a*b + c, could use FMA3 or FMA4 */
#undef gmx_nmsub_pr
#undef gmx_max_pr
#undef gmx_cmplt_pr
#undef gmx_and_pr
#undef gmx_or_pr
#undef gmx_andnot_pr

/* Only used for PBC in bonded interactions, can be avoided */
#undef gmx_round_pr
/* Not required, only used to speed up the nbnxn tabulated PME kernels */
#undef gmx_floor_pr
/* Not required, only used with x86 when blendv is faster than comparison */
#undef gmx_blendv_pr

#undef gmx_movemask_pr

/* Integer casts are only used for nbnxn x86 exclusion masks */
#undef gmx_mm_castsi128_pr
#undef gmx_mm_castsi256_pr

/* Conversions only used for nbnxn x86 exclusion masks and PME table lookup */
#undef gmx_cvttpr_epi32
#undef gmx_cvtepi32_pr

#undef gmx_invsqrt_pr
/* sqrt+inv+sin+cos+acos+atan2 are only used for bonded potentials */
#undef gmx_sqrt_pr
#undef gmx_inv_pr
#undef gmx_sincos_pr
#undef gmx_acos_pr
#undef gmx_atan_pr

#undef gmx_calc_rsq_pr
#undef gmx_sum4_pr

/* Only required for nbnxn analytical PME kernels */
#undef gmx_pmecorrF_pr
#undef gmx_pmecorrV_pr


/* Half SIMD-width types and operations only for nbnxn 2xnn search+kernels */
#undef gmx_mm_hpr

#undef gmx_load_hpr
#undef gmx_load1_hpr
#undef gmx_store_hpr
#undef gmx_add_hpr
#undef gmx_sub_hpr

#undef gmx_sum4_hpr

#undef gmx_2hpr_to_pr


/* By defining GMX_MM128_HERE or GMX_MM256_HERE before including this file
 * the same intrinsics, with defines, can be compiled for either 128 or 256
 * bit wide SIMD instructions.
 * On x86, the gmx_ prefix is replaced by _mm_ or _mm256_ (SSE or AVX-256).
 * On x86, the _pr suffix is replaced by _ps or _pd (single or double precision).
 * Note that on x86 compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */

#if !defined GMX_MM128_HERE && !defined GMX_MM256_HERE
#error "You should define GMX_MM128_HERE or GMX_MM256_HERE"
#endif

#if defined GMX_MM128_HERE && defined GMX_MM256_HERE
#error "You should not define both GMX_MM128_HERE and GMX_MM256_HERE"
#endif


#ifdef GMX_X86_SSE2

#ifdef GMX_MM128_HERE

#define gmx_epi32  __m128i

#ifndef GMX_DOUBLE

#include "gmx_x86_simd_single.h"

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m128

#define gmx_load_pr       _mm_load_ps
#define gmx_load1_pr      _mm_load1_ps
#define gmx_set1_pr       _mm_set1_ps
#define gmx_setzero_pr    _mm_setzero_ps
#define gmx_store_pr      _mm_store_ps
#define gmx_storeu_pr     _mm_storeu_ps

#define gmx_add_pr        _mm_add_ps
#define gmx_sub_pr        _mm_sub_ps
#define gmx_mul_pr        _mm_mul_ps
#ifdef GMX_X86_AVX_128_FMA
#define gmx_madd_pr(a, b, c)   _mm_macc_ps(a, b, c)
#define gmx_nmsub_pr(a, b, c)  _mm_nmacc_ps(a, b, c)
#else
#define gmx_madd_pr(a, b, c)   _mm_add_ps(c, _mm_mul_ps(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#endif
#define gmx_max_pr        _mm_max_ps
#define gmx_cmplt_pr      _mm_cmplt_ps
#define gmx_and_pr        _mm_and_ps
#define gmx_or_pr         _mm_or_ps
#define gmx_andnot_pr     _mm_andnot_ps

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_ps(x, 0x0)
#define gmx_floor_pr      _mm_floor_ps
#else
#define gmx_round_pr(x)   _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#endif
#if defined GMX_X86_AVX_128_FMA || defined GMX_X86_AVX_256
#define gmx_blendv_pr     _mm_blendv_ps
#endif

#define gmx_movemask_pr   _mm_movemask_ps

#define gmx_mm_castsi128_pr gmx_mm_castsi128_ps

#define gmx_cvttpr_epi32  _mm_cvttps_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_ps

#define gmx_invsqrt_pr    gmx_mm_invsqrt_ps
#define gmx_sqrt_pr       gmx_mm_sqrt_ps
#define gmx_inv_pr        gmx_mm_inv_ps
#define gmx_sincos_pr     gmx_mm_sincos_ps
#define gmx_acos_pr       gmx_mm_acos_ps
#define gmx_atan2_pr      gmx_mm_atan2_ps

#define gmx_calc_rsq_pr   gmx_mm_calc_rsq_ps
#define gmx_sum4_pr       gmx_mm_sum4_ps

#define gmx_pmecorrF_pr   gmx_mm_pmecorrF_ps
#define gmx_pmecorrV_pr   gmx_mm_pmecorrV_ps

#else /* ifndef GMX_DOUBLE */

#include "gmx_x86_simd_double.h"

#define GMX_SIMD_WIDTH_HERE  2

#define gmx_mm_pr  __m128d

#define gmx_load_pr       _mm_load_pd
#define gmx_load1_pr      _mm_load1_pd
#define gmx_set1_pr       _mm_set1_pd
#define gmx_setzero_pr    _mm_setzero_pd
#define gmx_store_pr      _mm_store_pd
#define gmx_storeu_pr     _mm_storeu_pd

#define gmx_add_pr        _mm_add_pd
#define gmx_sub_pr        _mm_sub_pd
#define gmx_mul_pr        _mm_mul_pd
#ifdef GMX_X86_AVX_128_FMA
#define gmx_madd_pr(a, b, c)   _mm_macc_pd(a, b, c)
#define gmx_nmsub_pr(a, b, c)  _mm_nmacc_pd(a, b, c)
#else
#define gmx_madd_pr(a, b, c)   _mm_add_pd(c, _mm_mul_pd(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm_sub_pd(c, _mm_mul_pd(a, b))
#endif
#define gmx_max_pr        _mm_max_pd
#define gmx_cmplt_pr      _mm_cmplt_pd
#define gmx_and_pr        _mm_and_pd
#define gmx_or_pr         _mm_or_pd
#define gmx_andnot_pr     _mm_andnot_pd

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_pd(x, 0x0)
#define gmx_floor_pr      _mm_floor_pd
#else
#define gmx_round_pr(x)   _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
#endif
#if defined GMX_X86_AVX_128_FMA || defined GMX_X86_AVX_256
#define gmx_blendv_pr     _mm_blendv_pd
#endif

#define gmx_movemask_pr   _mm_movemask_pd

#define gmx_mm_castsi128_pr gmx_mm_castsi128_pd

#define gmx_cvttpr_epi32  _mm_cvttpd_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_pd

#define gmx_invsqrt_pr    gmx_mm_invsqrt_pd
#define gmx_sqrt_pr       gmx_mm_sqrt_pd
#define gmx_inv_pr        gmx_mm_inv_pd
#define gmx_sincos_pr     gmx_mm_sincos_pd
#define gmx_acos_pr       gmx_mm_acos_pd
#define gmx_atan2_pr      gmx_mm_atan2_pd

#define gmx_calc_rsq_pr   gmx_mm_calc_rsq_pd
#define gmx_sum4_pr       gmx_mm_sum4_pd

#define gmx_pmecorrF_pr   gmx_mm_pmecorrF_pd
#define gmx_pmecorrV_pr   gmx_mm_pmecorrV_pd

#endif /* ifndef GMX_DOUBLE */

#endif /* GMX_MM128_HERE */

#ifdef GMX_MM256_HERE

#define gmx_epi32 __m256i

#ifndef GMX_DOUBLE

#include "gmx_x86_simd_single.h"

#define GMX_SIMD_WIDTH_HERE  8

#define gmx_mm_pr  __m256

#define gmx_load_pr       _mm256_load_ps
#define gmx_load1_pr(x)   _mm256_set1_ps((x)[0])
#define gmx_set1_pr       _mm256_set1_ps
#define gmx_setzero_pr    _mm256_setzero_ps
#define gmx_store_pr      _mm256_store_ps
#define gmx_storeu_pr     _mm256_storeu_ps

#define gmx_add_pr        _mm256_add_ps
#define gmx_sub_pr        _mm256_sub_ps
#define gmx_mul_pr        _mm256_mul_ps
#define gmx_madd_pr(a, b, c)   _mm256_add_ps(c, _mm256_mul_ps(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm256_sub_ps(c, _mm256_mul_ps(a, b))
#define gmx_max_pr        _mm256_max_ps
/* Not-equal (ordered, non-signaling)  */
#define gmx_cmpneq_pr(x, y)  _mm256_cmp_ps(x, y, 0x0c)
/* Less-than (ordered, non-signaling)  */
#define gmx_cmplt_pr(x, y) _mm256_cmp_ps(x, y, 0x11)
#define gmx_and_pr        _mm256_and_ps
#define gmx_or_pr         _mm256_or_ps
#define gmx_andnot_pr     _mm256_andnot_ps

#define gmx_round_pr(x)   _mm256_round_ps(x, 0x0)
#define gmx_floor_pr      _mm256_floor_ps
#define gmx_blendv_pr     _mm256_blendv_ps

#define gmx_movemask_pr   _mm256_movemask_ps

#define gmx_mm_castsi256_pr _mm256_castsi256_ps

#define gmx_cvttpr_epi32  _mm256_cvttps_epi32

#define gmx_invsqrt_pr    gmx_mm256_invsqrt_ps
#define gmx_sqrt_pr       gmx_mm256_sqrt_ps
#define gmx_inv_pr        gmx_mm256_inv_ps
#define gmx_sincos_pr     gmx_mm256_sincos_ps
#define gmx_acos_pr       gmx_mm256_acos_ps
#define gmx_atan2_pr      gmx_mm256_atan2_ps

#define gmx_calc_rsq_pr   gmx_mm256_calc_rsq_ps
#define gmx_sum4_pr       gmx_mm256_sum4_ps

#define gmx_pmecorrF_pr   gmx_mm256_pmecorrF_ps
#define gmx_pmecorrV_pr   gmx_mm256_pmecorrV_ps

#define gmx_loaddh_pr     gmx_mm256_load4_ps

/* Half SIMD-width type */
#define gmx_mm_hpr  __m128

/* Half SIMD-width macros */
#define gmx_load_hpr      _mm_load_ps
#define gmx_load1_hpr(x)  _mm_set1_ps((x)[0])
#define gmx_store_hpr     _mm_store_ps
#define gmx_add_hpr       _mm_add_ps
#define gmx_sub_hpr       _mm_sub_ps

#define gmx_sum4_hpr      gmx_mm256_sum4h_m128

/* Conversion between half and full SIMD-width */
#define gmx_2hpr_to_pr    gmx_mm256_set_m128

#else

#include "gmx_x86_simd_double.h"

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m256d

#define gmx_load_pr       _mm256_load_pd
#define gmx_load1_pr(x)   _mm256_set1_pd((x)[0])
#define gmx_set1_pr       _mm256_set1_pd
#define gmx_setzero_pr    _mm256_setzero_pd
#define gmx_store_pr      _mm256_store_pd
#define gmx_storeu_pr     _mm256_storeu_pd

#define gmx_add_pr        _mm256_add_pd
#define gmx_sub_pr        _mm256_sub_pd
#define gmx_mul_pr        _mm256_mul_pd
#define gmx_madd_pr(a, b, c)   _mm256_add_pd(c, _mm256_mul_pd(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#define gmx_max_pr        _mm256_max_pd
/* Not-equal (ordered, non-signaling)  */
#define gmx_cmpneq_pr(x, y)  _mm256_cmp_pd(x, y, 0x0c)
/* Less-than (ordered, non-signaling)  */
#define gmx_cmplt_pr(x, y) _mm256_cmp_pd(x, y, 0x11)
#define gmx_and_pr        _mm256_and_pd
#define gmx_or_pr         _mm256_or_pd
#define gmx_andnot_pr     _mm256_andnot_pd

#define gmx_round_pr(x)   _mm256_round_pd(x, 0x0)
#define gmx_floor_pr      _mm256_floor_pd
#define gmx_blendv_pr     _mm256_blendv_pd

#define gmx_movemask_pr   _mm256_movemask_pd

#define gmx_mm_castsi256_pr _mm256_castsi256_pd

#define gmx_cvttpr_epi32  _mm256_cvttpd_epi32

#define gmx_invsqrt_pr    gmx_mm256_invsqrt_pd
#define gmx_sqrt_pr       gmx_mm256_sqrt_pd
#define gmx_inv_pr        gmx_mm256_inv_pd
#define gmx_sincos_pr     gmx_mm256_sincos_pd
#define gmx_acos_pr       gmx_mm256_acos_pd
#define gmx_atan2_pr      gmx_mm256_atan2_pd

#define gmx_calc_rsq_pr   gmx_mm256_calc_rsq_pd
#define gmx_sum4_pr       gmx_mm256_sum4_pd

#define gmx_pmecorrF_pr   gmx_mm256_pmecorrF_pd
#define gmx_pmecorrV_pr   gmx_mm256_pmecorrV_pd

#endif

#endif /* GMX_MM256_HERE */

#endif /* GMX_X86_SSE2 */
