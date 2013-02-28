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

/* float/double SIMD register type */
#undef gmx_mm_pr

/* integer SIMD register type, only used in the tabulated PME kernels */
#undef gmx_epi32

#undef gmx_load_pr
#undef gmx_load1_pr
#undef gmx_set1_pr
#undef gmx_setzero_pr
#undef gmx_store_pr

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
/* gmx_blendzero_pr(real a, boolean b) does: (b ? a : 0) */
#undef gmx_blendzero_pr
/* Logical operations on SIMD booleans */
#undef gmx_and_pr
#undef gmx_or_pr
#undef gmx_andnot_pr

/* Only used for PBC in bonded interactions, can be avoided */
#undef gmx_round_pr
/* Not required, only used to speed up the nbnxn tabulated PME kernels */
#undef GMX_HAVE_SIMD_FLOOR
#undef gmx_floor_pr

/* Not required, only used when blendv is faster than comparison */
#undef GMX_HAVE_SIMD_BLENDV
#undef gmx_blendv_pr
/* Not required, gmx_anytrue(x) returns if any of the boolean is x is True.
 * If this is not present, define GMX_SIMD_IS_TRUE(real x),
 * which should return x==True, where True is True as defined in SIMD.
 */
#undef GMX_HAVE_SIMD_ANYTRUE
#undef gmx_anytrue_pr

/* Integer set and cast are only used for nbnxn exclusion masks */
#undef gmx_set1_epi32
#undef gmx_castsi_pr
/* For topology exclusion pair checking we need: ((a & b) ? True : False)
 * when we do a bit-wise and between a and b.
 * When integer SIMD operations are present, we use gmx_checkbitmask_epi32(a, b)
 * Otherwise we do all operations, except for the set1, in reals.
 */
#undef gmx_load_si
/* If the same bit is set in both input masks, return all bits 1, otherwise 0 */
#undef gmx_checkbitmask_epi32
/* As gmx_checkbitmask_epi32, but operates on reals. In double precision two
 * identical 32-bit masks are set in one double and one or both can be used.
 */
#undef gmx_checkbitmask_pr

/* Conversions only used for PME table lookup */
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


/* The same SIMD macros can be translated to SIMD intrinsics (and compiled
 * to instructions for) different SIMD width and float precision.
 *
 * On x86: The gmx_ prefix is replaced by _mm_ or _mm256_ (SSE or AVX).
 * The _pr suffix is replaced by _ps or _pd (for single or double precision).
 * Compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */


/* Generic macros for obtaining a SIMD aligned pointer from pointer x */
#undef gmx_simd_align_real
#undef gmx_simd_align_int


#ifdef GMX_USE_HALF_WIDTH_SIMD_HERE
#if defined GMX_X86_AVX_256
/* We have half SIMD width support, continue */
#else
#error "half SIMD width intrinsics are not supported"
#endif
#endif


#ifdef GMX_X86_SSE2

#if !defined GMX_X86_AVX_256 || defined GMX_USE_HALF_WIDTH_SIMD_HERE

#ifndef GMX_DOUBLE

#include "gmx_x86_simd_single.h"

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m128

#define gmx_epi32  __m128i

#define gmx_load_pr       _mm_load_ps
#define gmx_load1_pr      _mm_load1_ps
#define gmx_set1_pr       _mm_set1_ps
#define gmx_setzero_pr    _mm_setzero_ps
#define gmx_store_pr      _mm_store_ps

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
#define gmx_blendzero_pr  _mm_and_ps
#define gmx_and_pr        _mm_and_ps
#define gmx_or_pr         _mm_or_ps
#define gmx_andnot_pr     _mm_andnot_ps

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_ps(x, 0x0)
#define GMX_HAVE_SIMD_FLOOR
#define gmx_floor_pr      _mm_floor_ps
#else
#define gmx_round_pr(x)   _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#endif

#ifdef GMX_X86_SSE4_1
#define GMX_HAVE_SIMD_BLENDV
#define gmx_blendv_pr     _mm_blendv_ps
#endif

#define GMX_HAVE_SIMD_ANYTRUE
#define gmx_anytrue_pr    _mm_movemask_ps

#define gmx_set1_epi32    _mm_set1_epi32
#define gmx_castsi_pr     gmx_mm_castsi128_ps
#define gmx_load_si(i)    _mm_load_si128((__m128i *) (i))
#define gmx_checkbitmask_epi32(m0, m1) _mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128())

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

#define gmx_epi32  __m128i

#define gmx_load_pr       _mm_load_pd
#define gmx_load1_pr      _mm_load1_pd
#define gmx_set1_pr       _mm_set1_pd
#define gmx_setzero_pr    _mm_setzero_pd
#define gmx_store_pr      _mm_store_pd

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
#define gmx_blendzero_pr  _mm_and_pd
#define gmx_and_pr        _mm_and_pd
#define gmx_or_pr         _mm_or_pd
#define gmx_andnot_pr     _mm_andnot_pd

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_pd(x, 0x0)
#define GMX_HAVE_SIMD_FLOOR
#define gmx_floor_pr      _mm_floor_pd
#else
#define gmx_round_pr(x)   _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
/* gmx_floor_pr is not used in code for pre-SSE4_1 hardware */
#endif

#ifdef GMX_X86_SSE4_1
#define GMX_HAVE_SIMD_BLENDV
#define gmx_blendv_pr     _mm_blendv_pd
#endif

#define GMX_HAVE_SIMD_ANYTRUE
#define gmx_anytrue_pr    _mm_movemask_pd

#define gmx_set1_epi32    _mm_set1_epi32
#define gmx_castsi_pr     gmx_mm_castsi128_pd
#define gmx_load_si(i)    _mm_load_si128((__m128i *) (i))
#define gmx_checkbitmask_epi32(m0, m1) _mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128())

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

#else
/* We have GMX_X86_AVX_256 and not GMX_USE_HALF_WIDTH_SIMD_HERE,
 * so we use 256-bit SIMD.
 */

#ifndef GMX_DOUBLE

#include "gmx_x86_simd_single.h"

#define GMX_SIMD_WIDTH_HERE  8

#define gmx_mm_pr  __m256

#define gmx_epi32  __m256i

#define gmx_load_pr       _mm256_load_ps
#define gmx_load1_pr(x)   _mm256_set1_ps((x)[0])
#define gmx_set1_pr       _mm256_set1_ps
#define gmx_setzero_pr    _mm256_setzero_ps
#define gmx_store_pr      _mm256_store_ps

#define gmx_add_pr        _mm256_add_ps
#define gmx_sub_pr        _mm256_sub_ps
#define gmx_mul_pr        _mm256_mul_ps
#define gmx_madd_pr(a, b, c)   _mm256_add_ps(c, _mm256_mul_ps(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm256_sub_ps(c, _mm256_mul_ps(a, b))
#define gmx_max_pr        _mm256_max_ps
/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_cmplt_pr(x, y) _mm256_cmp_ps(x, y, 0x11)
#define gmx_blendzero_pr  _mm256_and_ps
#define gmx_and_pr        _mm256_and_ps
#define gmx_or_pr         _mm256_or_ps
#define gmx_andnot_pr     _mm256_andnot_ps

#define gmx_round_pr(x)   _mm256_round_ps(x, 0x0)
#define GMX_HAVE_SIMD_FLOOR
#define gmx_floor_pr      _mm256_floor_ps

#define GMX_HAVE_SIMD_BLENDV
#define gmx_blendv_pr     _mm256_blendv_ps

#define GMX_HAVE_SIMD_ANYTRUE
#define gmx_anytrue_pr    _mm256_movemask_ps

#define gmx_set1_epi32    _mm256_set1_epi32
#define gmx_castsi_pr     _mm256_castsi256_ps
/* With <= 16 bits used the cast and conversion should not be required,
 * since only mantissa bits are set and that would give a non-zero float,
 * but with the Intel compiler this does not work correctly.
 */
#define gmx_checkbitmask_pr(m0, m1) _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(m0, m1))), _mm256_setzero_ps(), 0x0c)

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

#else

#include "gmx_x86_simd_double.h"

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m256d

/* We use 128-bit integer registers because of missing 256-bit operations */
#define gmx_epi32  __m128i

#define gmx_load_pr       _mm256_load_pd
#define gmx_load1_pr(x)   _mm256_set1_pd((x)[0])
#define gmx_set1_pr       _mm256_set1_pd
#define gmx_setzero_pr    _mm256_setzero_pd
#define gmx_store_pr      _mm256_store_pd

#define gmx_add_pr        _mm256_add_pd
#define gmx_sub_pr        _mm256_sub_pd
#define gmx_mul_pr        _mm256_mul_pd
#define gmx_madd_pr(a, b, c)   _mm256_add_pd(c, _mm256_mul_pd(a, b))
#define gmx_nmsub_pr(a, b, c)  _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#define gmx_max_pr        _mm256_max_pd
/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_cmplt_pr(x, y) _mm256_cmp_pd(x, y, 0x11)
#define gmx_blendzero_pr  _mm256_and_pd
#define gmx_and_pr        _mm256_and_pd
#define gmx_or_pr         _mm256_or_pd
#define gmx_andnot_pr     _mm256_andnot_pd

#define gmx_round_pr(x)   _mm256_round_pd(x, 0x0)
#define GMX_HAVE_SIMD_FLOOR
#define gmx_floor_pr      _mm256_floor_pd

#define GMX_HAVE_SIMD_BLENDV
#define gmx_blendv_pr     _mm256_blendv_pd

#define GMX_HAVE_SIMD_ANYTRUE
#define gmx_anytrue_pr    _mm256_movemask_pd

#define gmx_set1_epi32    _mm256_set1_epi32
#define gmx_castsi_pr     _mm256_castsi256_pd
/* With <= 16 bits used the cast and conversion should not be required,
 * since only mantissa bits are set and that would give a non-zero float,
 * but with the Intel compiler this does not work correctly.
 * Because AVX does not have int->double conversion, we convert via float.
 */
#define gmx_checkbitmask_pr(m0, m1) _mm256_cmp_pd(_mm256_castps_pd(_mm256_cvtepi32_ps(_mm256_castpd_si256(_mm256_and_pd(m0, m1)))), _mm256_setzero_pd(), 0x0c)

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

#endif /* GMX_DOUBLE */

#endif /* 128- or 256-bit x86 SIMD */

#endif /* GMX_X86_SSE2 */


/* Generic macros to extract a SIMD aligned pointer from a pointer x.
 * x should have at least GMX_SIMD_WIDTH_HERE elements extra compared
 * to how many you want to use, to avoid indexing outside the aligned region.
 */

#define gmx_simd_align_real(x)  (real *)(((size_t)((x)+GMX_SIMD_WIDTH_HERE)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(real)-1))))

#define gmx_simd_align_int(x)   (int  *)(((size_t)((x)+GMX_SIMD_WIDTH_HERE)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(int )-1))))
