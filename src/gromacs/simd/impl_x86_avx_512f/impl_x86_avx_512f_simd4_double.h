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

#ifndef GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_DOUBLE_H

#include <math.h>

#include <immintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_avx_512f_common.h"
#include "impl_x86_avx_512f_simd_double.h"

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
/* We use __m256d to only access part of the registers, but for a few operations
 * we cast to the full register width when those operations are cheaper. We
 * also save some register space by using mask registers for the booleans.
 */
#define gmx_simd4_double_t          __m256d
#define gmx_simd4_load_d            _mm256_load_pd
#define gmx_simd4_load1_d           _mm256_broadcast_sd
#define gmx_simd4_set1_d            _mm256_set1_pd
#define gmx_simd4_store_d           _mm256_store_pd
#define gmx_simd4_loadu_d           _mm256_loadu_pd
#define gmx_simd4_storeu_d          _mm256_storeu_pd
#define gmx_simd4_setzero_d         _mm256_setzero_pd
#define gmx_simd4_add_d             _mm256_add_pd
#define gmx_simd4_sub_d             _mm256_sub_pd
#define gmx_simd4_mul_d             _mm256_mul_pd
#define gmx_simd4_fmadd_d           _mm256_fmadd_pd
#define gmx_simd4_fmsub_d           _mm256_fmsub_pd
#define gmx_simd4_fnmadd_d          _mm256_fnmadd_pd
#define gmx_simd4_fnmsub_d          _mm256_fnmsub_pd
#define gmx_simd4_and_d             _mm256_and_pd
#define gmx_simd4_andnot_d          _mm256_andnot_pd
#define gmx_simd4_or_d              _mm256_or_pd
#define gmx_simd4_xor_d             _mm256_xor_pd
/* We need to use the new table lookup instructions since we have specified
 * 14 bits of accuracy for AVX-512F.
 */
#define gmx_simd4_rsqrt_d(x)        _mm512_castpd512_pd256(_mm512_rsqrt14_pd(_mm512_castpd256_pd512(x)))
/* abs/neg cannot cause FP exceptions, so we can operate on entire register */
#define gmx_simd4_fabs_d(x)         _mm512_castpd512_pd256(_mm512_abs_pd(_mm512_castpd256_pd512(x)))
#define gmx_simd4_fneg_d(x)         _mm256_xor_pd(x, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
#define gmx_simd4_max_d(a, b)       _mm256_max_pd(a, b)
#define gmx_simd4_min_d(a, b)       _mm256_min_pd(a, b)
#define gmx_simd4_round_d(a)        _mm256_round_pd(a, _MM_FROUND_NINT)
#define gmx_simd4_trunc_d(a)        _mm256_round_pd(a, _MM_FROUND_TRUNC)
#define gmx_simd4_dotproduct3_d(a, b) gmx_simd4_dotproduct3_d_x86_avx_512f(a, b)
#define gmx_simd4_dbool_t           __mmask16
#define gmx_simd4_cmpeq_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_EQ_OQ)
#define gmx_simd4_cmplt_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_LT_OS)
#define gmx_simd4_cmple_d(a, b)     _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b), _CMP_LE_OS)
#define gmx_simd4_and_db            _mm512_kand
#define gmx_simd4_or_db             _mm512_kor
#define gmx_simd4_anytrue_db(x)     (_mm512_mask2int(x)&0xF)
#define gmx_simd4_blendzero_d(a, sel)    _mm512_castpd512_pd256(_mm512_mask_mov_pd(_mm512_setzero_pd(), sel, _mm512_castpd256_pd512(a)))
#define gmx_simd4_blendnotzero_d(a, sel) _mm512_castpd512_pd256(_mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(sel), _mm512_castpd256_pd512(a)))
#define gmx_simd4_blendv_d(a, b, sel)    _mm512_castpd512_pd256(_mm512_mask_blend_pd(sel, _mm512_castpd256_pd512(a), _mm512_castpd256_pd512(b)))
#define gmx_simd4_reduce_d(x)       gmx_simd4_reduce_d_x86_avx_512f(x)

/* Implementation helpers */

static gmx_inline double
gmx_simd4_reduce_d_x86_avx_512f(__m256d a)
{
    double  f;
    __m128d a0, a1;
    a  = _mm256_hadd_pd(a, a);
    a0 = _mm256_castpd256_pd128(a);
    a1 = _mm256_extractf128_pd(a, 0x1);
    a0 = _mm_add_sd(a0, a1);
    _mm_store_sd(&f, a0);
    return f;
}

static gmx_inline double
gmx_simd4_dotproduct3_d_x86_avx_512f(__m256d a, __m256d b)
{
    double  d;
    __m128d tmp1, tmp2;
    a    = _mm256_mul_pd(a, b);
    tmp1 = _mm256_castpd256_pd128(a);
    tmp2 = _mm256_extractf128_pd(a, 0x1);

    tmp1 = _mm_add_pd(tmp1, _mm_permute_pd(tmp1, _MM_SHUFFLE2(0, 1)));
    tmp1 = _mm_add_pd(tmp1, tmp2);
    _mm_store_sd(&d, tmp1);
    return d;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512F_SIMD4_DOUBLE_H */
