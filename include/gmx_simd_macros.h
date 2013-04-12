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

#ifdef _gmx_simd_macros_h_
#error "gmx_simd_macros.h included twice"
#else
#define _gmx_simd_macros_h_

/* NOTE: SSE2 acceleration does not include floor or blendv */


/* Uncomment the next line, without other SIMD active, for testing plain-C */
/* #define GMX_SIMD_PLAIN_C */
#ifdef GMX_SIMD_PLAIN_C
/* Plain C 2-way SIMD reference implementation, also serves as documentation */
#define GMX_HAVE_SIMD_MACROS

#define GMX_SIMD_WIDTH_HERE  2

/* float/double SIMD register type */
typedef struct {
    real x;
    real y;
} gmx_mm_pr;

/* boolean SIMD register type */
typedef struct {
    char x;
    char y;
} gmx_mm_pb;

/* integer SIMD register type, only used in tabulated non-bonded kernels */
typedef struct {
    unsigned x;
    unsigned y;
} gmx_epi32;
#define GMX_SIMD_EPI32_WIDTH  2

/* Load GMX_SIMD_WIDTH_HERE reals for memory starting at r */
static gmx_inline gmx_mm_pr gmx_load_pr(const real *r) { gmx_mm_pr a; a.x = r[0]; a.y = r[1]; return a; };
/* Set all SIMD register elements to *r */
static gmx_inline gmx_mm_pr gmx_load1_pr(const real *r) { gmx_mm_pr a; a.x = *r; a.y = *r; return a; };
static gmx_inline gmx_mm_pr gmx_set1_pr(real r) { gmx_mm_pr a; a.x = r; a.y = r; return a; };
static gmx_inline gmx_mm_pr gmx_setzero_pr() { gmx_mm_pr a; a.x = 0.0; a.y = 0.0; return a; };
static gmx_inline void gmx_store_pr(real *dest, gmx_mm_pr src) { dest[0] = src.x; dest[1] = src.y; };

static gmx_inline gmx_mm_pr gmx_add_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = a.x + b.x; c.y = a.y + b.y; return c; };
static gmx_inline gmx_mm_pr gmx_sub_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = a.x - b.x; c.y = a.y - b.y; return c; };
static gmx_inline gmx_mm_pr gmx_mul_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = a.x*b.x; c.y = a.y*b.y; return c; };
/* For the FMA macros below, aim for c=d in code, so FMA3 uses 1 instruction */
static gmx_inline gmx_mm_pr gmx_madd_pr(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c) { gmx_mm_pr d; d.x = a.x*b.x + c.x; d.y = a.y*b.y + c.y; return d; };
static gmx_inline gmx_mm_pr gmx_nmsub_pr(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c) { gmx_mm_pr d; d.x = -a.x*b.x + c.x; d.y = -a.y*b.y + c.y; return d; };
static gmx_inline gmx_mm_pr gmx_max_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = (a.x >= b.x) ? a.x : b.x; c.y = (a.y >= b.y) ? a.y : b.y; return c; };
static gmx_inline gmx_mm_pr gmx_blendzero_pr(gmx_mm_pr a, gmx_mm_pb b) { gmx_mm_pr c; c.x = b.x ? a.x : 0.0; c.y = b.y ? a.y : 0.0; return c; };

static gmx_inline gmx_mm_pr gmx_round_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = (real)(int)(a.x + (a.x >= 0 ? 0.5 : -0.5)); b.y = (real)(int)(a.y + (a.y >= 0 ? 0.5 : -0.5)); return b; };
/* Not required, only used to speed up the nbnxn tabulated PME kernels */
#define GMX_SIMD_HAVE_FLOOR
#ifdef GMX_SIMD_HAVE_FLOOR
static gmx_inline gmx_mm_pr gmx_floor_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = (real)(int)((a.x >= 0) ? a.x : a.x - 1); b.y = (real)(int)((a.y >= 0) ? a.y : a.y - 1); return b; };
#endif

/* Not required, only used when blendv is faster than comparison */
#define GMX_SIMD_HAVE_BLENDV
#ifdef GMX_SIMD_HAVE_BLENDV
static gmx_inline gmx_mm_pr gmx_blendv_pr(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c) { gmx_mm_pr d; d.x = (c.x >= 0) ? a.x : b.x; d.y = (c.y >= 0) ? a.y : b.y; return d; };
#endif

/* Copy the sign of a to b, assumes b >= 0 for efficiency */
static gmx_inline gmx_mm_pr gmx_cpsgn_nonneg_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = (a.x >= 0) ? b.x : -b.x; c.y = (a.y >= 0) ? b.y : -b.y; return c; };

/* Very specific operation required in the non-bonded kernels */
static gmx_inline gmx_mm_pr gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_pr b, gmx_mm_pr c) { gmx_mm_pr d; d.x = a.x ? b.x : b.x + c.x; d.y = a.y ? b.y : b.y + c.y; return d; };

/* Comparison */
static gmx_inline gmx_mm_pb gmx_cmplt_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pb c; c.x = (a.x < b.x); c.y = (a.y < b.y); return c; };

/* Logical operations on SIMD booleans */
static gmx_inline gmx_mm_pb gmx_and_pb(gmx_mm_pb a, gmx_mm_pb b) { gmx_mm_pb c; c.x = a.x && b.x; c.y = a.y && b.y; return c; };
static gmx_inline gmx_mm_pb gmx_or_pb(gmx_mm_pb a, gmx_mm_pb b) { gmx_mm_pb c; c.x = a.x || b.x; c.y = a.y || b.y; return c; };

/* Not required, gmx_anytrue(x) returns if any of the boolean is x is True.
 * If this is not present, define GMX_SIMD_IS_TRUE(real x),
 * which should return x==True, where True is True as defined in SIMD.
 */
#define GMX_SIMD_HAVE_ANYTRUE
#ifdef GMX_SIMD_HAVE_ANYTRUE
static gmx_inline int gmx_anytrue_pb(gmx_mm_pb a) { return (a.x != 0) || (a.y != 0); };
#else
/* If we don't have gmx_anytrue_pb, we need to store gmx_mm_pb */
static gmx_inline void gmx_store_pb(real *dest, gmx_mm_pb src) { dest[0] = src.x; dest[1] = src.y; };
#endif

/* For topology exclusion pair checking we need: ((a & b) ? True : False)
 * when we do a bit-wise and between a and b.
 * When integer SIMD operations are present, we use gmx_checkbitmask_epi32(a, b)
 * Otherwise we do all operations, except for the set1, in reals.
 */

#define GMX_SIMD_HAVE_CHECKBITMASK_EPI32
#ifdef GMX_SIMD_HAVE_CHECKBITMASK_EPI32

/* Integer set and cast are only used for nbnxn exclusion masks */
static gmx_inline gmx_epi32 gmx_set1_epi32(unsigned i) { gmx_epi32 a; a.x = i; a.y = i; return a; };

static gmx_inline gmx_epi32 gmx_load_si(unsigned *i) { gmx_epi32 a; a.x = i[0]; a.y = i[1]; return a; };
/* If the same bit is set in both input masks, return TRUE, else FALSE */

static gmx_inline gmx_mm_pb gmx_checkbitmask_epi32(gmx_epi32 a, gmx_epi32 b) { gmx_mm_pb c; c.x = (a.x & b.x) != 0; c.y = (a.y & b.y) != 0; return c; };

#endif

/* #define GMX_SIMD_HAVE_CHECKBITMASK_PR */
#ifdef GMX_SIMD_HAVE_CHECKBITMASK_PR
static gmx_inline gmx_mm_pr gmx_castsi_pr(gmx_epi32 a) { };
/* As gmx_checkbitmask_epi32, but operates on reals. In double precision two
 * identical 32-bit masks are set in one double and one or both can be used.
 */
static gmx_inline gmx_mm_pb gmx_checkbitmask_pr(gmx_mm_pr a, gmx_mm_pr b) { };
#endif

/* Conversions only used for PME table lookup */
static gmx_inline gmx_epi32 gmx_cvttpr_epi32(gmx_mm_pr a) { gmx_epi32 b; b.x = (unsigned)a.x; b.y = (unsigned)a.y; return b; };
static gmx_inline gmx_mm_pr gmx_cvtepi32_pr(gmx_epi32 a) { gmx_mm_pr b; b.x = (real)a.x; b.y = (real)a.y; return b; };

/* These two function only need to be approximate, Newton-Raphson iteration
 * is used for full accuracy in gmx_invsqrt_pr and gmx_inv_pr.
 */
static gmx_inline gmx_mm_pr gmx_rsqrt_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = 1.0/sqrt(a.x); b.y = 1.0/sqrt(a.y); return b; };
static gmx_inline gmx_mm_pr gmx_rcp_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = 1.0/a.x; b.y = 1.0/a.y; return b; };

/* sqrt+inv+sin+cos+acos+atan2 are used for bonded potentials, exp for PME */
#define GMX_SIMD_HAVE_EXP
#ifdef GMX_SIMD_HAVE_EXP
static gmx_inline gmx_mm_pr gmx_exp_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = exp(a.x); b.y = exp(a.y); return b; };
#endif
#define GMX_SIMD_HAVE_TRIGONOMETRIC
#ifdef GMX_SIMD_HAVE_TRIGONOMETRIC
static gmx_inline gmx_mm_pr gmx_sqrt_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = sqrt(a.x); b.y = sqrt(a.y); return b; };
static gmx_inline void gmx_sincos_pr(gmx_mm_pr a, gmx_mm_pr *s, gmx_mm_pr *c) { s->x = sin(a.x); s->y = sin(a.y); c->x = cos(a.x); c->y = cos(a.y); };
static gmx_inline gmx_mm_pr gmx_acos_pr(gmx_mm_pr a) { gmx_mm_pr b; b.x = acos(a.x); b.y = acos(a.y); return b; };
static gmx_inline gmx_mm_pr gmx_atan2_pr(gmx_mm_pr a, gmx_mm_pr b) { gmx_mm_pr c; c.x = atan2(a.x, b.x); c.y = atan2(a.y, b.y); return c; };
#endif

#endif /* GMX_SIMD_PLAIN_C */


/* The same SIMD macros can be translated to SIMD intrinsics (and compiled
 * to instructions for) different SIMD width and float precision.
 *
 * On x86: The gmx_ prefix is replaced by _mm_ or _mm256_ (SSE or AVX).
 * The _pr suffix is replaced by _ps or _pd (for single or double precision).
 * Compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */


#ifdef GMX_USE_HALF_WIDTH_SIMD_HERE
#if defined GMX_X86_AVX_256
/* We have half SIMD width support, continue */
#else
#error "half SIMD width intrinsics are not supported"
#endif
#endif


#ifdef GMX_X86_SSE2
/* This is for general x86 SIMD instruction sets that also support SSE2 */
#define GMX_HAVE_SIMD_MACROS

/* Include the highest supported x86 SIMD intrisics + math functions */
#ifdef GMX_X86_AVX_256
#include "gmx_x86_avx_256.h"
#ifdef GMX_DOUBLE
#include "gmx_math_x86_avx_256_double.h"
#else
#include "gmx_math_x86_avx_256_single.h"
#endif
#else
#ifdef GMX_X86_AVX_128_FMA
#include "gmx_x86_avx_128_fma.h"
#ifdef GMX_DOUBLE
#include "gmx_math_x86_avx_128_fma_double.h"
#else
#include "gmx_math_x86_avx_128_fma_single.h"
#endif
#else
#ifdef GMX_X86_SSE4_1
#include "gmx_x86_sse4_1.h"
#ifdef GMX_DOUBLE
#include "gmx_math_x86_sse4_1_double.h"
#else
#include "gmx_math_x86_sse4_1_single.h"
#endif
#else
#ifdef GMX_X86_SSE2
#include "gmx_x86_sse2.h"
#ifdef GMX_DOUBLE
#include "gmx_math_x86_sse2_double.h"
#else
#include "gmx_math_x86_sse2_single.h"
#endif
#else
#error No x86 acceleration defined
#endif
#endif
#endif
#endif
/* exp and trigonometric functions are included above */
#define GMX_SIMD_HAVE_EXP
#define GMX_SIMD_HAVE_TRIGONOMETRIC

#if !defined GMX_X86_AVX_256 || defined GMX_USE_HALF_WIDTH_SIMD_HERE

#ifndef GMX_DOUBLE

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m128

#define gmx_mm_pb  __m128

#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

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
#define gmx_blendzero_pr  _mm_and_ps

#define gmx_cmplt_pr      _mm_cmplt_ps
#define gmx_and_pb        _mm_and_ps
#define gmx_or_pb         _mm_or_ps

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_ps(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_floor_pr      _mm_floor_ps
#else
#define gmx_round_pr(x)   _mm_cvtepi32_ps(_mm_cvtps_epi32(x))
#endif

#ifdef GMX_X86_SSE4_1
#define GMX_SIMD_HAVE_BLENDV
#define gmx_blendv_pr     _mm_blendv_ps
#endif

static gmx_inline gmx_mm_pr gmx_cpsgn_nonneg_pr(gmx_mm_pr a, gmx_mm_pr b)
{
    /* The value -0.0 has only the sign-bit set */
    gmx_mm_pr sign_mask = _mm_set1_ps(-0.0);
    return _mm_or_ps(_mm_and_ps(a, sign_mask), b);
};

static gmx_inline gmx_mm_pr gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_pr b, gmx_mm_pr c) { return _mm_add_ps(b, _mm_andnot_ps(a, c)); };

#define GMX_SIMD_HAVE_ANYTRUE
#define gmx_anytrue_pb    _mm_movemask_ps

#define GMX_SIMD_HAVE_CHECKBITMASK_EPI32
#define gmx_set1_epi32    _mm_set1_epi32
#define gmx_load_si(i)    _mm_load_si128((__m128i *) (i))
#define gmx_checkbitmask_epi32(m0, m1) gmx_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()))

#define gmx_cvttpr_epi32  _mm_cvttps_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_ps

#define gmx_rsqrt_pr      _mm_rsqrt_ps
#define gmx_rcp_pr        _mm_rcp_ps

#define gmx_exp_pr        gmx_mm_exp_ps
#define gmx_sqrt_pr       gmx_mm_sqrt_ps
#define gmx_sincos_pr     gmx_mm_sincos_ps
#define gmx_acos_pr       gmx_mm_acos_ps
#define gmx_atan2_pr      gmx_mm_atan2_ps

#else /* ifndef GMX_DOUBLE */

#define GMX_SIMD_WIDTH_HERE  2

#define gmx_mm_pr  __m128d

#define gmx_mm_pb  __m128d

#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

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
#define gmx_blendzero_pr  _mm_and_pd

#ifdef GMX_X86_SSE4_1
#define gmx_round_pr(x)   _mm_round_pd(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_floor_pr      _mm_floor_pd
#else
#define gmx_round_pr(x)   _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
/* gmx_floor_pr is not used in code for pre-SSE4_1 hardware */
#endif

#ifdef GMX_X86_SSE4_1
#define GMX_SIMD_HAVE_BLENDV
#define gmx_blendv_pr     _mm_blendv_pd
#endif

static gmx_inline gmx_mm_pr gmx_cpsgn_nonneg_pr(gmx_mm_pr a, gmx_mm_pr b)
{
    gmx_mm_pr sign_mask = _mm_set1_pd(-0.0);
    return _mm_or_pd(_mm_and_pd(a, sign_mask), b);
};

static gmx_inline gmx_mm_pr gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_pr b, gmx_mm_pr c) { return _mm_add_pd(b, _mm_andnot_pd(a, c)); };

#define gmx_cmplt_pr      _mm_cmplt_pd

#define gmx_and_pb        _mm_and_pd
#define gmx_or_pb         _mm_or_pd

#define GMX_SIMD_HAVE_ANYTRUE
#define gmx_anytrue_pb    _mm_movemask_pd

#define GMX_SIMD_HAVE_CHECKBITMASK_EPI32
#define gmx_set1_epi32    _mm_set1_epi32
#define gmx_load_si(i)    _mm_load_si128((__m128i *) (i))
#define gmx_checkbitmask_epi32(m0, m1) gmx_mm_castsi128_pd(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()))

#define gmx_cvttpr_epi32  _mm_cvttpd_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_pd

#define gmx_rsqrt_pr(r)   _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(r)))
#define gmx_rcp_pr(r)     _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(r)))

#define gmx_exp_pr        gmx_mm_exp_pd
#define gmx_sqrt_pr       gmx_mm_sqrt_pd
#define gmx_sincos_pr     gmx_mm_sincos_pd
#define gmx_acos_pr       gmx_mm_acos_pd
#define gmx_atan2_pr      gmx_mm_atan2_pd

#endif /* ifndef GMX_DOUBLE */

#else
/* We have GMX_X86_AVX_256 and not GMX_USE_HALF_WIDTH_SIMD_HERE,
 * so we use 256-bit SIMD.
 */

#ifndef GMX_DOUBLE

#define GMX_SIMD_WIDTH_HERE  8

#define gmx_mm_pr  __m256

#define gmx_mm_pb  __m256

#define gmx_epi32  __m256i
#define GMX_SIMD_EPI32_WIDTH  8

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
#define gmx_blendzero_pr  _mm256_and_ps

#define gmx_round_pr(x)   _mm256_round_ps(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_floor_pr      _mm256_floor_ps

#define GMX_SIMD_HAVE_BLENDV
#define gmx_blendv_pr     _mm256_blendv_ps

static gmx_inline gmx_mm_pr gmx_cpsgn_nonneg_pr(gmx_mm_pr a, gmx_mm_pr b)
{
    gmx_mm_pr sign_mask = _mm256_set1_ps(-0.0);
    return _mm256_or_ps(_mm256_and_ps(a, sign_mask), b);
};

static gmx_inline gmx_mm_pr gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_pr b, gmx_mm_pr c) { return _mm256_add_ps(b, _mm256_andnot_ps(a, c)); };

/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_cmplt_pr(x, y) _mm256_cmp_ps(x, y, 0x11)
#define gmx_and_pb        _mm256_and_ps
#define gmx_or_pb         _mm256_or_ps

#define GMX_SIMD_HAVE_ANYTRUE
#define gmx_anytrue_pb    _mm256_movemask_ps

#define GMX_SIMD_HAVE_CHECKBITMASK_PR
#define gmx_set1_epi32    _mm256_set1_epi32
#define gmx_castsi_pr     _mm256_castsi256_ps
/* With <= 16 bits used the cast and conversion should not be required,
 * since only mantissa bits are set and that would give a non-zero float,
 * but with the Intel compiler this does not work correctly.
 */
#define gmx_checkbitmask_pr(m0, m1) _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(m0, m1))), _mm256_setzero_ps(), 0x0c)

#define gmx_cvttpr_epi32  _mm256_cvttps_epi32

#define gmx_rsqrt_pr      _mm256_rsqrt_ps
#define gmx_rcp_pr        _mm256_rcp_ps

#define gmx_exp_pr        gmx_mm256_exp_ps
#define gmx_sqrt_pr       gmx_mm256_sqrt_ps
#define gmx_sincos_pr     gmx_mm256_sincos_ps
#define gmx_acos_pr       gmx_mm256_acos_ps
#define gmx_atan2_pr      gmx_mm256_atan2_ps

#else

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m256d

#define gmx_mm_pb  __m256d

/* We use 128-bit integer registers because of missing 256-bit operations */
#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

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
#define gmx_blendzero_pr  _mm256_and_pd

#define gmx_round_pr(x)   _mm256_round_pd(x, 0x0)
#define GMX_SIMD_HAVE_FLOOR
#define gmx_floor_pr      _mm256_floor_pd

#define GMX_SIMD_HAVE_BLENDV
#define gmx_blendv_pr     _mm256_blendv_pd

static gmx_inline gmx_mm_pr gmx_cpsgn_nonneg_pr(gmx_mm_pr a, gmx_mm_pr b)
{
    gmx_mm_pr sign_mask = _mm256_set1_pd(-0.0);
    return _mm256_or_pd(_mm256_and_pd(a, sign_mask), b);
};

static gmx_inline gmx_mm_pr gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_pr b, gmx_mm_pr c) { return _mm256_add_pd(b, _mm256_andnot_pd(a, c)); };

/* Less-than (we use ordered, non-signaling, but that's not required) */
#define gmx_cmplt_pr(x, y) _mm256_cmp_pd(x, y, 0x11)

#define gmx_and_pb        _mm256_and_pd
#define gmx_or_pb         _mm256_or_pd

#define GMX_SIMD_HAVE_ANYTRUE
#define gmx_anytrue_pb    _mm256_movemask_pd

#define GMX_SIMD_HAVE_CHECKBITMASK_PR
#define gmx_set1_epi32    _mm256_set1_epi32
#define gmx_castsi_pr     _mm256_castsi256_pd
/* With <= 16 bits used the cast and conversion should not be required,
 * since only mantissa bits are set and that would give a non-zero float,
 * but with the Intel compiler this does not work correctly.
 * Because AVX does not have int->double conversion, we convert via float.
 */
#define gmx_checkbitmask_pr(m0, m1) _mm256_cmp_pd(_mm256_castps_pd(_mm256_cvtepi32_ps(_mm256_castpd_si256(_mm256_and_pd(m0, m1)))), _mm256_setzero_pd(), 0x0c)

#define gmx_cvttpr_epi32  _mm256_cvttpd_epi32

#define gmx_rsqrt_pr(r)   _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(r)))
#define gmx_rcp_pr(r)     _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(r)))

#define gmx_exp_pr        gmx_mm256_exp_pd
#define gmx_sqrt_pr       gmx_mm256_sqrt_pd
#define gmx_sincos_pr     gmx_mm256_sincos_pd
#define gmx_acos_pr       gmx_mm256_acos_pd
#define gmx_atan2_pr      gmx_mm256_atan2_pd

#endif /* GMX_DOUBLE */

#endif /* 128- or 256-bit x86 SIMD */

#endif /* GMX_X86_SSE2 */


#ifdef GMX_HAVE_SIMD_MACROS
/* Generic functions to extract a SIMD aligned pointer from a pointer x.
 * x should have at least GMX_SIMD_WIDTH_HERE elements extra compared
 * to how many you want to use, to avoid indexing outside the aligned region.
 */

static gmx_inline real *
gmx_simd_align_real(const real *x)
{
    return (real *)(((size_t)((x)+GMX_SIMD_WIDTH_HERE)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(real)-1))));
}

static gmx_inline int *
gmx_simd_align_int(const int *x)
{
    return (int  *)(((size_t)((x)+GMX_SIMD_WIDTH_HERE)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(int )-1))));
}


/* Include the math functions which only need the above macros,
 * generally these are the ones that don't need masking operations.
 */
#ifdef GMX_DOUBLE
#include "gmx_simd_math_double.h"
#else
#include "gmx_simd_math_single.h"
#endif

#endif /* GMX_HAVE_SIMD_MACROS */

#endif /* _gmx_simd_macros_h_ */
