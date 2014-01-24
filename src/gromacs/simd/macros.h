/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

/* The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding macros here should be (nearly)
 * all that is needed.
 */

#ifdef GMX_SIMD_MACROS_H
#error "gromacs/simd/macros.h included twice"
#else
#define GMX_SIMD_MACROS_H

/* NOTE: SSE2 acceleration does not include floor or blendv */


/* Uncomment the next line, without other SIMD active, for testing plain-C */
/* #define GMX_SIMD_REFERENCE */
#ifdef GMX_SIMD_REFERENCE
/* Plain C SIMD reference implementation, also serves as documentation */
#define GMX_HAVE_SIMD_MACROS

/* In general the reference SIMD supports any SIMD width, including 1.
 * See types/nb_verlet.h for details
 */
#define GMX_SIMD_REF_WIDTH  4

/* Include plain-C reference implementation, also serves as documentation */
#include "gromacs/simd/macros_ref.h"

#define GMX_SIMD_REAL_WIDTH  GMX_SIMD_REF_WIDTH

/* float/double SIMD register type */
#define gmx_simd_real_t  gmx_simd_ref_pr

/* boolean SIMD register type */
#define gmx_simd_bool_t  gmx_simd_ref_pb

/* integer SIMD register type, only for table indexing and exclusion masks */
#define gmx_simd_int32_t  gmx_simd_ref_epi32
#define GMX_SIMD_INT32_WIDTH  GMX_SIMD_REF_EPI32_WIDTH

/* Load GMX_SIMD_REAL_WIDTH reals for memory starting at r */
#define gmx_simd_load_r       gmx_simd_ref_load_pr
/* Set all SIMD register elements to *r */
#define gmx_simd_load1_r      gmx_simd_ref_load1_pr
#define gmx_simd_set1_r       gmx_simd_ref_set1_pr
#define gmx_simd_setzero_r    gmx_simd_ref_setzero_pr
#define gmx_simd_store_r      gmx_simd_ref_store_pr

#define gmx_simd_add_r        gmx_simd_ref_add_pr
#define gmx_simd_sub_r        gmx_simd_ref_sub_pr
#define gmx_simd_mul_r        gmx_simd_ref_mul_pr
/* For the FMA macros below, aim for c=d in code, so FMA3 uses 1 instruction */
#define gmx_simd_fmadd_r       gmx_simd_ref_madd_pr
#define gmx_simd_fnmadd_r      gmx_simd_ref_nmsub_pr

#define gmx_simd_max_r        gmx_simd_ref_max_pr
#define gmx_simd_blendzero_r  gmx_simd_ref_blendzero_pr

#define gmx_simd_round_r      gmx_simd_ref_round_pr

/* Not required, only used to speed up the nbnxn tabulated PME kernels */
#define GMX_SIMD_HAVE_FLOOR
#ifdef GMX_SIMD_HAVE_FLOOR
#define gmx_simd_floor_r      gmx_simd_ref_floor_pr
#endif

/* Not required, only used when blendv is faster than comparison */
#define GMX_SIMD_HAVE_BLENDV
#ifdef GMX_SIMD_HAVE_BLENDV
#define gmx_simd_blendv_r     gmx_simd_ref_blendv_pr
#endif

/* Copy the sign of a to b, assumes b >= 0 for efficiency */
#define gmx_cpsgn_nonneg_pr  gmx_simd_ref_cpsgn_nonneg_pr

/* Very specific operation required in the non-bonded kernels */
#define gmx_masknot_add_pr   gmx_simd_ref_masknot_add_pr

/* Comparison */
#define gmx_simd_cmplt_r      gmx_simd_ref_cmplt_pr

/* Logical operations on SIMD booleans */
#define gmx_simd_and_b        gmx_simd_ref_and_pb
#define gmx_simd_or_b         gmx_simd_ref_or_pb

/* Returns a single int (0/1) which tells if any of the 4 booleans is True */
#define gmx_simd_anytrue_b    gmx_simd_ref_anytrue_pb

/* Conversions only used for PME table lookup */
#define gmx_simd_cvtt_r2i  gmx_simd_ref_cvttpr_epi32
#define gmx_simd_cvt_i2r   gmx_simd_ref_cvtepi32_pr

/* These two function only need to be approximate, Newton-Raphson iteration
 * is used for full accuracy in gmx_simd_invsqrt_r and gmx_simd_inv_r.
 */
#define gmx_simd_rsqrt_r      gmx_simd_ref_rsqrt_pr
#define gmx_simd_rcp_r        gmx_simd_ref_rcp_pr

/* sqrt+inv+sin+cos+acos+atan2 are used for bonded potentials, exp for PME */
#define GMX_SIMD_HAVE_EXP
#ifdef GMX_SIMD_HAVE_EXP
#define gmx_simd_exp_r        gmx_simd_ref_exp_pr
#endif
#define GMX_SIMD_HAVE_TRIGONOMETRIC
#ifdef GMX_SIMD_HAVE_TRIGONOMETRIC
#define gmx_simd_sqrt_r       gmx_simd_ref_sqrt_pr
#define gmx_simd_sincos_r     gmx_simd_ref_sincos_pr
#define gmx_simd_acos_r       gmx_simd_ref_acos_pr
#define gmx_simd_atan2_r      gmx_simd_ref_atan2_pr
#endif

#endif /* GMX_SIMD_REFERENCE */


/* The same SIMD macros can be translated to SIMD intrinsics (and compiled
 * to instructions for) different SIMD width and float precision.
 *
 * On x86: The gmx_ prefix is replaced by _mm_ or _mm256_ (SSE or AVX).
 * The _pr suffix is replaced by _ps or _pd (for single or double precision).
 * Compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */


#ifdef GMX_USE_HALF_WIDTH_SIMD_HERE
#if defined GMX_X86_AVX_256_OR_HIGHER || defined __MIC__
/* We have half SIMD width support, continue */
#else
#error "half SIMD width intrinsics are not supported"
#endif
#endif

#if defined GMX_TARGET_X86 && !defined __MIC__

#ifdef GMX_X86_SSE2_OR_HIGHER
/* This is for general x86 SIMD instruction sets that also support SSE2 */
#define GMX_HAVE_SIMD_MACROS

/* Include the highest supported x86 SIMD intrisics + math functions */
#ifdef GMX_X86_AVX_256_OR_HIGHER
#include "general_x86_avx_256.h"
#ifdef GMX_DOUBLE
#include "math_x86_avx_256_double.h"
#else  /* GMX_DOUBLE */
#include "math_x86_avx_256_single.h"
#endif /* GMX_DOUBLE */
#else  /* GMX_X86_AVX_256_OR_HIGHER */
#ifdef GMX_X86_AVX_128_FMA_OR_HIGHER
#include "general_x86_avx_128_fma.h"
#ifdef GMX_DOUBLE
#include "math_x86_avx_128_fma_double.h"
#else  /* GMX_DOUBLE */
#include "math_x86_avx_128_fma_single.h"
#endif /* GMX_DOUBLE */
#else  /* GMX_X86_AVX_128_FMA_OR_HIGHER */
#ifdef GMX_SIMD_X86_SSE4_1
#include "general_x86_sse4_1.h"
#ifdef GMX_DOUBLE
#include "math_x86_sse4_1_double.h"
#else  /* GMX_DOUBLE */
#include "math_x86_sse4_1_single.h"
#endif /* GMX_DOUBLE */
#else  /* GMX_X86_SSE4_1_OR_HIGHER */
#ifdef GMX_X86_SSE2_OR_HIGHER
#include "general_x86_sse2.h"
#ifdef GMX_DOUBLE
#include "math_x86_sse2_double.h"
#else  /* GMX_DOUBLE */
#include "math_x86_sse2_single.h"
#endif /* GMX_DOUBLE */
#else  /* GMX_X86_SSE2_OR_HIGHER */
#error No x86 acceleration defined
#endif /* GMX_X86_SSE2_OR_HIGHER */
#endif /* GMX_X86_SSE4_1_OR_HIGHER */
#endif /* GMX_X86_AVX_128_FMA_OR_HIGHER */
#endif /* GMX_X86_AVX_256_OR_HIGHER */

/* exp and trigonometric functions are included above */
#define GMX_SIMD_HAVE_EXP
#define GMX_SIMD_HAVE_TRIGONOMETRIC

#if !defined GMX_X86_AVX_256_OR_HIGHER || defined GMX_USE_HALF_WIDTH_SIMD_HERE

#ifndef GMX_DOUBLE

#define GMX_SIMD_REAL_WIDTH  4

#define gmx_simd_real_t  __m128

#define gmx_simd_bool_t  __m128

#define gmx_simd_int32_t  __m128i
#define GMX_SIMD_INT32_WIDTH  4

#define gmx_simd_load_r       _mm_load_ps
#define gmx_simd_load1_r      _mm_load1_ps
#define gmx_simd_set1_r       _mm_set1_ps
#define gmx_simd_setzero_r    _mm_setzero_ps
#define gmx_simd_store_r      _mm_store_ps

#define gmx_simd_add_r        _mm_add_ps
#define gmx_simd_sub_r        _mm_sub_ps
#define gmx_simd_mul_r        _mm_mul_ps
#ifdef GMX_X86_AVX_128_FMA_OR_HIGHER
#define GMX_SIMD_HAVE_FMA
#define gmx_simd_fmadd_r(a, b, c)   _mm_macc_ps(a, b, c)
#define gmx_simd_fnmadd_r(a, b, c)  _mm_nmacc_ps(a, b, c)
#else
#define gmx_simd_fmadd_r(a, b, c)   _mm_add_ps(c, _mm_mul_ps(a, b))
#define gmx_simd_fnmadd_r(a, b, c)  _mm_sub_ps(c, _mm_mul_ps(a, b))
#endif
#define gmx_simd_max_r        _mm_max_ps
#define gmx_simd_blendzero_r  _mm_and_ps

#define gmx_simd_cmplt_r      _mm_cmplt_ps
#define gmx_simd_and_b        _mm_and_ps
#define gmx_simd_or_b         _mm_or_ps

#ifdef GMX_X86_SSE4_1_OR_HIGHER
#define gmx_simd_round_r(x)   _mm_round_ps(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_simd_floor_r      _mm_floor_ps
#else
#define gmx_simd_round_r(x)   _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#endif

#ifdef GMX_X86_SSE4_1_OR_HIGHER
#define GMX_SIMD_HAVE_BLENDV
#define gmx_simd_blendv_r     _mm_blendv_ps
#endif

static gmx_inline gmx_simd_real_t gmx_cpsgn_nonneg_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    /* The value -0.0 has only the sign-bit set */
    gmx_simd_real_t sign_mask = _mm_set1_ps(-0.0);
    return _mm_or_ps(_mm_and_ps(a, sign_mask), b);
};

static gmx_inline gmx_simd_real_t gmx_masknot_add_pr(gmx_simd_bool_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return _mm_add_ps(b, _mm_andnot_ps(a, c));
};

#define gmx_simd_anytrue_b    _mm_movemask_ps

#define gmx_simd_cvtt_r2i  _mm_cvttps_epi32
#define gmx_simd_cvt_i2r   _mm_cvtepi32_ps

#define gmx_simd_rsqrt_r      _mm_rsqrt_ps
#define gmx_simd_rcp_r        _mm_rcp_ps

#define gmx_simd_exp_r        gmx_mm_exp_ps
#define gmx_simd_sqrt_r       gmx_mm_sqrt_ps
#define gmx_simd_sincos_r     gmx_mm_sincos_ps
#define gmx_simd_acos_r       gmx_mm_acos_ps
#define gmx_simd_atan2_r      gmx_mm_atan2_ps
#define gmx_simd_erfc_r       gmx_mm_erfc_ps

#else /* ifndef GMX_DOUBLE */

#define GMX_SIMD_REAL_WIDTH  2

#define gmx_simd_real_t  __m128d

#define gmx_simd_bool_t  __m128d

#define gmx_simd_int32_t  __m128i
#define GMX_SIMD_INT32_WIDTH  4

#define gmx_simd_load_r       _mm_load_pd
#define gmx_simd_load1_r      _mm_load1_pd
#define gmx_simd_set1_r       _mm_set1_pd
#define gmx_simd_setzero_r    _mm_setzero_pd
#define gmx_simd_store_r      _mm_store_pd

#define gmx_simd_add_r        _mm_add_pd
#define gmx_simd_sub_r        _mm_sub_pd
#define gmx_simd_mul_r        _mm_mul_pd
#ifdef GMX_X86_AVX_128_FMA_OR_HIGHER
#define GMX_SIMD_HAVE_FMA
#define gmx_simd_fmadd_r(a, b, c)   _mm_macc_pd(a, b, c)
#define gmx_simd_fnmadd_r(a, b, c)  _mm_nmacc_pd(a, b, c)
#else
#define gmx_simd_fmadd_r(a, b, c)   _mm_add_pd(c, _mm_mul_pd(a, b))
#define gmx_simd_fnmadd_r(a, b, c)  _mm_sub_pd(c, _mm_mul_pd(a, b))
#endif
#define gmx_simd_max_r        _mm_max_pd
#define gmx_simd_blendzero_r  _mm_and_pd

#ifdef GMX_X86_SSE4_1_OR_HIGHER
#define gmx_simd_round_r(x)   _mm_round_pd(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_simd_floor_r      _mm_floor_pd
#else
#define gmx_simd_round_r(x)   _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
/* gmx_simd_floor_r is not used in code for pre-SSE4_1 hardware */
#endif

#ifdef GMX_X86_SSE4_1_OR_HIGHER
#define GMX_SIMD_HAVE_BLENDV
#define gmx_simd_blendv_r     _mm_blendv_pd
#endif

static gmx_inline gmx_simd_real_t gmx_cpsgn_nonneg_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    gmx_simd_real_t sign_mask = _mm_set1_pd(-0.0);
    return _mm_or_pd(_mm_and_pd(a, sign_mask), b);
};

static gmx_inline gmx_simd_real_t gmx_masknot_add_pr(gmx_simd_bool_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return _mm_add_pd(b, _mm_andnot_pd(a, c));
};

#define gmx_simd_cmplt_r      _mm_cmplt_pd

#define gmx_simd_and_b        _mm_and_pd
#define gmx_simd_or_b         _mm_or_pd

#define gmx_simd_anytrue_b    _mm_movemask_pd

#define gmx_simd_cvtt_r2i  _mm_cvttpd_epi32
#define gmx_simd_cvt_i2r   _mm_cvtepi32_pd

#define gmx_simd_rsqrt_r(r)   _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(r)))
#define gmx_simd_rcp_r(r)     _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(r)))

#define gmx_simd_exp_r        gmx_mm_exp_pd
#define gmx_simd_sqrt_r       gmx_mm_sqrt_pd
#define gmx_simd_sincos_r     gmx_mm_sincos_pd
#define gmx_simd_acos_r       gmx_mm_acos_pd
#define gmx_simd_atan2_r      gmx_mm_atan2_pd
#define gmx_simd_erfc_r       gmx_mm_erfc_pd

#endif /* ifndef GMX_DOUBLE */

#else
/* We have GMX_X86_AVX_256_OR_HIGHER and not GMX_USE_HALF_WIDTH_SIMD_HERE,
 * so we use 256-bit SIMD.
 */

#ifndef GMX_DOUBLE

#define GMX_SIMD_REAL_WIDTH  8

#define gmx_simd_real_t  __m256

#define gmx_simd_bool_t  __m256

#define gmx_simd_int32_t  __m256i
#define GMX_SIMD_INT32_WIDTH  8

#define gmx_simd_load_r       _mm256_load_ps
#define gmx_simd_load1_r(x)   _mm256_set1_ps((x)[0])
#define gmx_simd_set1_r       _mm256_set1_ps
#define gmx_simd_setzero_r    _mm256_setzero_ps
#define gmx_simd_store_r      _mm256_store_ps

#define gmx_simd_add_r        _mm256_add_ps
#define gmx_simd_sub_r        _mm256_sub_ps
#define gmx_simd_mul_r        _mm256_mul_ps
#define gmx_simd_fmadd_r(a, b, c)   _mm256_add_ps(c, _mm256_mul_ps(a, b))
#define gmx_simd_fnmadd_r(a, b, c)  _mm256_sub_ps(c, _mm256_mul_ps(a, b))
#define gmx_simd_max_r        _mm256_max_ps
#define gmx_simd_blendzero_r  _mm256_and_ps

#define gmx_simd_round_r(x)   _mm256_round_ps(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_simd_floor_r      _mm256_floor_ps

#define GMX_SIMD_HAVE_BLENDV
#define gmx_simd_blendv_r     _mm256_blendv_ps

static gmx_inline gmx_simd_real_t gmx_cpsgn_nonneg_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    gmx_simd_real_t sign_mask = _mm256_set1_ps(-0.0);
    return _mm256_or_ps(_mm256_and_ps(a, sign_mask), b);
};

static gmx_inline gmx_simd_real_t gmx_masknot_add_pr(gmx_simd_bool_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return _mm256_add_ps(b, _mm256_andnot_ps(a, c));
};

/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_simd_cmplt_r(x, y) _mm256_cmp_ps(x, y, 0x11)
#define gmx_simd_and_b        _mm256_and_ps
#define gmx_simd_or_b         _mm256_or_ps

#define gmx_simd_anytrue_b    _mm256_movemask_ps

#define gmx_simd_cvtt_r2i  _mm256_cvttps_epi32

#define gmx_simd_rsqrt_r      _mm256_rsqrt_ps
#define gmx_simd_rcp_r        _mm256_rcp_ps

#define gmx_simd_exp_r        gmx_mm256_exp_ps
#define gmx_simd_sqrt_r       gmx_mm256_sqrt_ps
#define gmx_simd_sincos_r     gmx_mm256_sincos_ps
#define gmx_simd_acos_r       gmx_mm256_acos_ps
#define gmx_simd_atan2_r      gmx_mm256_atan2_ps
#define gmx_simd_erfc_r       gmx_mm256_erfc_ps

#else /* ifndef GMX_DOUBLE */

#define GMX_SIMD_REAL_WIDTH  4

#define gmx_simd_real_t  __m256d

#define gmx_simd_bool_t  __m256d

/* We use 128-bit integer registers because of missing 256-bit operations */
#define gmx_simd_int32_t  __m128i
#define GMX_SIMD_INT32_WIDTH  4

#define gmx_simd_load_r       _mm256_load_pd
#define gmx_simd_load1_r(x)   _mm256_set1_pd((x)[0])
#define gmx_simd_set1_r       _mm256_set1_pd
#define gmx_simd_setzero_r    _mm256_setzero_pd
#define gmx_simd_store_r      _mm256_store_pd

#define gmx_simd_add_r        _mm256_add_pd
#define gmx_simd_sub_r        _mm256_sub_pd
#define gmx_simd_mul_r        _mm256_mul_pd
#define gmx_simd_fmadd_r(a, b, c)   _mm256_add_pd(c, _mm256_mul_pd(a, b))
#define gmx_simd_fnmadd_r(a, b, c)  _mm256_sub_pd(c, _mm256_mul_pd(a, b))
#define gmx_simd_max_r        _mm256_max_pd
#define gmx_simd_blendzero_r  _mm256_and_pd

#define gmx_simd_round_r(x)   _mm256_round_pd(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_simd_floor_r      _mm256_floor_pd

#define GMX_SIMD_HAVE_BLENDV
#define gmx_simd_blendv_r     _mm256_blendv_pd

static gmx_inline gmx_simd_real_t gmx_cpsgn_nonneg_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    gmx_simd_real_t sign_mask = _mm256_set1_pd(-0.0);
    return _mm256_or_pd(_mm256_and_pd(a, sign_mask), b);
};

static gmx_inline gmx_simd_real_t gmx_masknot_add_pr(gmx_simd_bool_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return _mm256_add_pd(b, _mm256_andnot_pd(a, c));
};

/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_simd_cmplt_r(x, y) _mm256_cmp_pd(x, y, 0x11)

#define gmx_simd_and_b        _mm256_and_pd
#define gmx_simd_or_b         _mm256_or_pd

#define gmx_simd_anytrue_b    _mm256_movemask_pd

#define gmx_simd_cvtt_r2i  _mm256_cvttpd_epi32

#define gmx_simd_rsqrt_r(r)   _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(r)))
#define gmx_simd_rcp_r(r)     _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(r)))

#define gmx_simd_exp_r        gmx_mm256_exp_pd
#define gmx_simd_sqrt_r       gmx_mm256_sqrt_pd
#define gmx_simd_sincos_r     gmx_mm256_sincos_pd
#define gmx_simd_acos_r       gmx_mm256_acos_pd
#define gmx_simd_atan2_r      gmx_mm256_atan2_pd
#define gmx_simd_erfc_r       gmx_mm256_erfc_pd

#endif /* ifndef GMX_DOUBLE */

#endif /* 128- or 256-bit x86 SIMD */

#endif /* GMX_X86_SSE2_OR_HIGHER */

#endif /* GMX_TARGET_X86 */

#ifdef GMX_SIMD_IBM_QPX

/* This hack works on the compilers that can reach this code. A real
   solution with broader scope will be proposed in master branch. */
#define gmx_always_inline __attribute__((always_inline))

/* This is for the A2 core on BlueGene/Q that supports IBM's QPX
   vector built-in functions */
#include <mass_simd.h>
#define GMX_HAVE_SIMD_MACROS
#ifdef __clang__
#include <qpxmath.h>
#endif

/* No need to version the code by the precision, because the QPX AXU
   extends to and truncates from double precision for free. */

#define GMX_SIMD_REAL_WIDTH  4
typedef vector4double gmx_simd_real_t;
typedef vector4double gmx_simd_bool_t;
typedef vector4double gmx_simd_int32_t;
#define GMX_SIMD_INT32_WIDTH  4

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_load_r(const real *a)
{
#ifdef NDEBUG
    return vec_ld(0, (real *) a);
#else
    return vec_lda(0, (real *) a);
#endif
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_load1_r(const real *a)
{
    return vec_splats(*a);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_set1_r(real a)
{
    return vec_splats(a);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_setzero_r()
{
    return vec_splats(0.0);
}

static gmx_inline void gmx_always_inline gmx_simd_store_r(real *a, gmx_simd_real_t b)
{
#ifdef NDEBUG
    vec_st(b, 0, a);
#else
    vec_sta(b, 0, a);
#endif
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_add_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_add(a, b);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_sub_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_sub(a, b);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_mul_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_mul(a, b);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_fmadd_r(gmx_simd_real_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return vec_madd(a, b, c);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_fnmadd_r(gmx_simd_real_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return vec_nmsub(a, b, c);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_max_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_sel(b, a, vec_sub(a, b));
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_blendzero_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_sel(gmx_simd_setzero_r(), a, b);
}

static gmx_inline gmx_simd_bool_t gmx_always_inline gmx_simd_cmplt_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_cmplt(a, b);
}

static gmx_inline gmx_simd_bool_t gmx_always_inline gmx_simd_and_b(gmx_simd_bool_t a, gmx_simd_bool_t b)
{
    return vec_and(a, b);
}

static gmx_inline gmx_simd_bool_t gmx_always_inline gmx_simd_or_b(gmx_simd_bool_t a, gmx_simd_bool_t b)
{
    return vec_or(a, b);
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_round_r(gmx_simd_real_t a)
{
    return vec_round(a);
}

#define GMX_SIMD_HAVE_FLOOR
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_floor_r(gmx_simd_real_t a)
{
    return vec_floor(a);
}

#define GMX_SIMD_HAVE_BLENDV
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_blendv_r(gmx_simd_real_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return vec_sel(b, a, gmx_simd_cmplt_r(gmx_simd_setzero_r(), c));
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_cpsgn_nonneg_pr(gmx_simd_real_t a, gmx_simd_real_t b)
{
    return vec_cpsgn(a, b);
};

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_masknot_add_pr(gmx_simd_bool_t a, gmx_simd_real_t b, gmx_simd_real_t c)
{
    return vec_add(b, vec_sel(c, gmx_simd_setzero_r(), a));
};

static gmx_inline gmx_bool gmx_always_inline
GMX_SIMD_IS_TRUE(real x)
{
    return x >= 0.0;
}

static gmx_inline gmx_simd_int32_t gmx_always_inline gmx_simd_cvtt_r2i(gmx_simd_real_t a)
{
    return vec_ctiwuz(a);
}
/* Don't want this, we have floor */
/* #define gmx_simd_cvt_i2r   vec_cvtepi32 */

/* A2 core on BG/Q delivers relative error of 2^-14, whereas Power ISA
   Architecture only promises 2^-8. So probably no need for
   Newton-Raphson iterates at single or double. */
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_rsqrt_r(gmx_simd_real_t a)
{
    return vec_rsqrte(a);
}

/* A2 core on BG/Q delivers relative error of 2^-14, whereas Power ISA
   Architecture only promises 2^-5. So probably no need for
   Newton-Raphson iterates at single or double. */
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_rcp_r(gmx_simd_real_t a)
{
    return vec_re(a);
}

/* Note that here, and below, we use the built-in SLEEF port when
   compiling on BlueGene/Q with clang */

#define GMX_SIMD_HAVE_EXP
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_exp_r(gmx_simd_real_t a)
{
#ifdef __clang__
#ifndef GMX_DOUBLE
    return xexpf(a);
#else
    return xexp(a);
#endif
#else
#ifndef GMX_DOUBLE
    return expf4(a);
#else
    return expd4(a);
#endif
#endif
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_sqrt_r(gmx_simd_real_t a)
{
#ifdef NDEBUG
    return vec_swsqrt_nochk(a);
#else
    return vec_swsqrt(a);
#endif
}

#define GMX_SIMD_HAVE_TRIGONOMETRIC
static gmx_inline int gmx_always_inline gmx_simd_sincos_r(gmx_simd_real_t a, gmx_simd_real_t *b, gmx_simd_real_t *c)
{
#ifdef __clang__
#ifndef GMX_DOUBLE
    xsincosf(a, b, c);
#else
    xsincos(a, b, c);
#endif
#else
#ifndef GMX_DOUBLE
    sincosf4(a, b, c);
#else
    sincosd4(a, b, c);
#endif
#endif
    return 1;
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_acos_r(gmx_simd_real_t a)
{
#ifdef __clang__
#ifndef GMX_DOUBLE
    return xacosf(a);
#else
    return xacos(a);
#endif
#else
#ifndef GMX_DOUBLE
    return acosf4(a);
#else
    return acosd4(a);
#endif
#endif
}

/* NB The order of parameters here is correct; the
   documentation of atan2[df]4 in SIMD MASS is wrong. */
static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_atan2_r(gmx_simd_real_t a, gmx_simd_real_t b)
{
#ifdef __clang__
#ifndef GMX_DOUBLE
    return xatan2f(a, b);
#else
    return xatan2(a, b);
#endif
#else
#ifndef GMX_DOUBLE
    return atan2f4(a, b);
#else
    return atan2d4(a, b);
#endif
#endif
}

static gmx_inline gmx_simd_real_t gmx_always_inline gmx_simd_erfc_r(gmx_simd_real_t a)
{
    /* The BG/Q qpxmath.h vector math library intended for use with
       bgclang does not have erfc, so we need to use a function from
       mass_simd.h. If this changes, then the #include <mass_simd.h> can
       become conditional. */
#ifndef GMX_DOUBLE
    return erfcf4(a);
#else
    return erfcd4(a);
#endif
}

/* TODO: gmx_mm_erfc_p[sd] should be generalized using gmx_*_pr, so that it just works on BlueGene */

static gmx_inline int gmx_always_inline
gmx_simd_anytrue_b(gmx_simd_bool_t a)
{
    /* The "anytrue" is done solely on the QPX AXU (which is the only
       available FPU). This is awkward, because pretty much no
       "horizontal" SIMD-vector operations exist, unlike x86 where
       SSE4.1 added various kinds of horizontal operations. So we have
       to make do with shifting vector elements and operating on the
       results. This makes for lots of data dependency, but the main
       alternative of storing to memory and reloading is not going to
       help, either. OpenMP over 2 or 4 hardware threads per core will
       hide much of the latency from the data dependency. The
       vec_extract() lets the compiler correctly use a floating-point
       comparison on the zeroth vector element, which avoids needing
       memory at all.
     */
    gmx_simd_bool_t vec_shifted_left_0 = a;
    gmx_simd_bool_t vec_shifted_left_1 = vec_sldw(a, a, 1);
    gmx_simd_bool_t vec_shifted_left_2 = vec_sldw(a, a, 2);
    gmx_simd_bool_t vec_shifted_left_3 = vec_sldw(a, a, 3);

    gmx_simd_bool_t vec_return = vec_or(vec_or(vec_shifted_left_2, vec_shifted_left_3),
                                        vec_or(vec_shifted_left_0, vec_shifted_left_1));
    return (0.0 < vec_extract(vec_return, 0));
};

#undef gmx_always_inline

#endif /* GMX_SIMD_IBM_QPX */

#ifdef __MIC__
#include "general_x86_mic.h"
#endif

#ifdef GMX_HAVE_SIMD_MACROS
/* Generic functions to extract a SIMD aligned pointer from a pointer x.
 * x should have at least GMX_SIMD_REAL_WIDTH elements extra compared
 * to how many you want to use, to avoid indexing outside the aligned region.
 */

static gmx_inline real *
gmx_simd_align_r(const real *x)
{
    return (real *)(((size_t)((x)+GMX_SIMD_REAL_WIDTH)) & (~((size_t)(GMX_SIMD_REAL_WIDTH*sizeof(real)-1))));
}

static gmx_inline int *
gmx_simd_align_int(const int *x)
{
    return (int  *)(((size_t)((x)+GMX_SIMD_REAL_WIDTH)) & (~((size_t)(GMX_SIMD_REAL_WIDTH*sizeof(int )-1))));
}


/* Include the math functions which only need the above macros,
 * generally these are the ones that don't need masking operations.
 */
#ifdef GMX_DOUBLE
#include "math_double.h"
#else
#include "math_single.h"
#endif


#endif /* GMX_HAVE_SIMD_MACROS */

#endif
