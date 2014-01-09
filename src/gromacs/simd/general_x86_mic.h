/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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

#ifndef _general_x86_mic_h_
#define _general_x86_mic_h_

/* This file contains the SIMD implmenetation for Intel MIC
 */

#include <math.h>
#include <immintrin.h>

#ifdef GMX_DOUBLE
#error "Double precision isn't supported on Intel Phi yet"
#endif

typedef __m512 gmx_mm_ps;
typedef __m512 gmx_mm_pr;
/* boolean SIMD register type */
typedef __mmask16 gmx_mm_pb;
typedef __m512i gmx_epi32;

#define GMX_HAVE_SIMD_MACROS
#define GMX_SIMD_WIDTH_HERE  16
#define GMX_SIMD_EPI32_WIDTH 16

#define gmx_load_pr _mm512_load_ps

/* Set all SIMD register elements to *r */
static gmx_inline gmx_mm_ps
gmx_load1_pr(const real *r)
{
    return _mm512_extload_ps(r, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
}

#define gmx_set1_pr _mm512_set1_ps
/* Set all SIMD register elements to 0 */
#define gmx_setzero_pr _mm512_setzero_ps
#define gmx_store_pr _mm512_store_ps

#define gmx_add_pr _mm512_add_ps
#define gmx_sub_pr _mm512_sub_ps
#define gmx_mul_pr _mm512_mul_ps

#define GMX_SIMD_HAVE_FMA
#define gmx_madd_pr _mm512_fmadd_ps
#define gmx_nmsub_pr _mm512_fnmadd_ps

#define gmx_max_pr _mm512_max_ps

static gmx_inline gmx_mm_ps
gmx_blendzero_pr(gmx_mm_ps a, gmx_mm_pb b)
{
    return _mm512_mask_mov_ps(_mm512_setzero_ps(), b, a);
}

#define gmx_round_pr _mm512_rint_ps

#define GMX_SIMD_HAVE_FLOOR
#define gmx_floor_pr _mm512_floor_ps

/* Copy the sign of a to b, assumes b >= 0 for efficiency */
static gmx_inline gmx_mm_ps
gmx_cpsgn_nonneg_pr(gmx_mm_ps a, gmx_mm_ps b)
{
    __m512 zero = _mm512_setzero_ps();
    __m512 neg1 = _mm512_set1_ps(-1);
    /* TODO (only bond): Bitwise operations on floating points can be done after casting to int.
       That allows us to do it the same way as AVX which might be faster. */
    return _mm512_mask_mul_ps(b, _mm512_cmplt_ps_mask(a, zero), b, neg1);
}

/* Very specific operation required in the non-bonded kernels */
static gmx_inline gmx_mm_ps
gmx_masknot_add_pr(gmx_mm_pb a, gmx_mm_ps b, gmx_mm_ps c)
{
    return _mm512_mask_add_ps(b, _mm512_knot(a), b, c);
}

/* Comparison */
#define gmx_cmplt_pr _mm512_cmplt_ps_mask

/* Logical AND on SIMD booleans. */
#define gmx_and_pb _mm512_kand

/* Logical OR on SIMD booleans. */
#define gmx_or_pb _mm512_kor

/* Returns a single int (0/1) which tells if any of the booleans is True
   It returns the full mask (not 1 for True). But given that any non-zero is True this is OK. */
#define gmx_anytrue_pb _mm512_mask2int

/* Conversions only used for PME table lookup */
static gmx_inline gmx_epi32
gmx_cvttpr_epi32(gmx_mm_ps a)
{
    return _mm512_cvtfxpnt_round_adjustps_epi32(a, _MM_ROUND_MODE_DOWN, _MM_EXPADJ_NONE);
};

/* These two function only need to be approximate, Newton-Raphson iteration
 * is used for full accuracy in gmx_invsqrt_pr and gmx_inv_pr.
 */
#define gmx_rsqrt_pr _mm512_rsqrt23_ps
#define gmx_rcp_pr _mm512_rcp23_ps

#define GMX_SIMD_HAVE_EXP
#define gmx_exp_pr _mm512_exp_ps

#define GMX_SIMD_HAVE_ERFC
#define gmx_erfc_pr _mm512_erfc_ps

#define GMX_SIMD_HAVE_TRIGONOMETRIC
#define gmx_sqrt_pr  _mm512_sqrt_ps

static gmx_inline int
gmx_sincos_pr(gmx_mm_ps a,
              gmx_mm_ps *s, gmx_mm_ps *c)
{
    /* TODO (only bond): optimize that both are calculated together.
       Or (if if that isn't fast on MIC) don't call sincos if only one is needed. */
    *s = _mm512_sin_ps(a);
    *c = _mm512_cos_ps(a);
    return 0;
}

#define gmx_acos_pr _mm512_acos_ps
#define gmx_atan2_pr _mm512_atan2_ps

#endif /* _general_x86_mic_h_ */
