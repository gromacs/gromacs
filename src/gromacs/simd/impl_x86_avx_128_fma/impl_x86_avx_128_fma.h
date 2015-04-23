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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_H

#include "config.h"

#include <assert.h>
#include <math.h>

#include <immintrin.h>
#include <x86intrin.h>

/* x86 128-bit AVX with FMA SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for details
 */

/* Inherit parts of AVX_128_FMA from SSE4.1 */
#include "gromacs/simd/impl_x86_sse4_1/impl_x86_sse4_1.h"
/* Increment over SSE4.1 capabilities */
#define GMX_SIMD_X86_AVX_128_FMA_OR_HIGHER

/* Override some capability definitions for things added in AVX over SSE4.1 */
#define GMX_SIMD_HAVE_FMA
#define GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD4_HAVE_DOUBLE  /* We can use 256-bit operations for this */


/* Work around gcc bug with wrong type for mask formal parameter to maskload/maskstore */
#ifdef GMX_SIMD_X86_AVX_GCC_MASKLOAD_BUG
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), _mm_castsi128_ps(mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)   _mm_maskstore_ps((mem), _mm_castsi128_ps(mask), (x))
#else
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), (mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)   _mm_maskstore_ps((mem), (mask), (x))
#endif


/* SINGLE */
#undef  gmx_simd_fmadd_f
#define gmx_simd_fmadd_f                 _mm_macc_ps
#undef  gmx_simd_fmsub_f
#define gmx_simd_fmsub_f                 _mm_msub_ps
#undef  gmx_simd_fnmadd_f
#define gmx_simd_fnmadd_f                _mm_nmacc_ps
#undef  gmx_simd_fnmsub_f
#define gmx_simd_fnmsub_f                _mm_nmsub_ps
#undef  gmx_simd_fraction_f
#define gmx_simd_fraction_f              _mm_frcz_ps
#undef  gmx_simd_gather_loadu_transpose_f
#define gmx_simd_gather_loadu_transpose_f       gmx_simd_gather_loadu_transpose_f_avx_128_fma
#undef  gmx_simd_transpose_scatter_storeu_f
#define gmx_simd_transpose_scatter_storeu_f     gmx_simd_transpose_scatter_storeu_f_avx_128_fma
#undef  gmx_simd_expand_scalars_to_triplets_f
#define gmx_simd_expand_scalars_to_triplets_f   gmx_simd_expand_scalars_to_triplets_f_avx_128_fma

/* DOUBLE */
#undef  gmx_simd_fmadd_d
#define gmx_simd_fmadd_d                _mm_macc_pd
#undef  gmx_simd_fmsub_d
#define gmx_simd_fmsub_d                 _mm_msub_pd
#undef  gmx_simd_fnmadd_d
#define gmx_simd_fnmadd_d                _mm_nmacc_pd
#undef  gmx_simd_fnmsub_d
#define gmx_simd_fnmsub_d                _mm_nmsub_pd
#undef  gmx_simd_fraction_d
#define gmx_simd_fraction_d              _mm_frcz_pd
#undef  gmx_simd_expand_scalars_to_triplets_d
#define gmx_simd_expand_scalars_to_triplets_d   gmx_simd_expand_scalars_to_triplets_d_avx_128_fma

/* Even if the _main_ SIMD implementation for this architecture file corresponds
 * to 128-bit AVX (since it will be faster), the 256-bit operations will always
 * be available in AVX, so we can use them for double precision SIMD4!
 */
/* SIMD4 Double precision floating point */
#define gmx_simd4_double_t               __m256d
#define gmx_simd4_load_d                 _mm256_load_pd
#define gmx_simd4_load1_d                _mm256_broadcast_sd
#define gmx_simd4_set1_d                 _mm256_set1_pd
#define gmx_simd4_store_d                _mm256_store_pd
#define gmx_simd4_loadu_d                _mm256_loadu_pd
#define gmx_simd4_storeu_d               _mm256_storeu_pd
#define gmx_simd4_setzero_d              _mm256_setzero_pd
#define gmx_simd4_add_d                  _mm256_add_pd
#define gmx_simd4_sub_d                  _mm256_sub_pd
#define gmx_simd4_mul_d                  _mm256_mul_pd
#define gmx_simd4_fmadd_d                _mm256_macc_pd
#define gmx_simd4_fmsub_d                _mm256_msub_pd
#define gmx_simd4_fnmadd_d               _mm256_nmacc_pd
#define gmx_simd4_fnmsub_d               _mm256_nmsub_pd
#define gmx_simd4_and_d                  _mm256_and_pd
#define gmx_simd4_andnot_d               _mm256_andnot_pd
#define gmx_simd4_or_d                   _mm256_or_pd
#define gmx_simd4_xor_d                  _mm256_xor_pd
#define gmx_simd4_rsqrt_d(x)             _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)))
#define gmx_simd4_fabs_d(x)              _mm256_andnot_pd(_mm256_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define gmx_simd4_fneg_d(x)              _mm256_xor_pd(x, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
#define gmx_simd4_max_d                  _mm256_max_pd
#define gmx_simd4_min_d                  _mm256_min_pd
#define gmx_simd4_round_d(x)             _mm256_round_pd(x, _MM_FROUND_NINT)
#define gmx_simd4_trunc_d(x)             _mm256_round_pd(x, _MM_FROUND_TRUNC)
#define gmx_simd4_dotproduct3_d          gmx_simd4_dotproduct3_d_avx_128_fma
#define gmx_simd4_transpose_d            gmx_simd4_transpose_d_avx_128_fma
/* SIMD4 booleans corresponding to double */
#define gmx_simd4_dbool_t                __m256d
#define gmx_simd4_cmpeq_d(a, b)           _mm256_cmp_pd(a, b, _CMP_EQ_OQ)
#define gmx_simd4_cmplt_d(a, b)           _mm256_cmp_pd(a, b, _CMP_LT_OQ)
#define gmx_simd4_cmple_d(a, b)           _mm256_cmp_pd(a, b, _CMP_LE_OQ)
#define gmx_simd4_and_db                 _mm256_and_pd
#define gmx_simd4_or_db                  _mm256_or_pd
#define gmx_simd4_anytrue_db             _mm256_movemask_pd
#define gmx_simd4_blendzero_d            _mm256_and_pd
#define gmx_simd4_blendnotzero_d(a, sel)  _mm256_andnot_pd(sel, a)
#define gmx_simd4_blendv_d               _mm256_blendv_pd
#define gmx_simd4_reduce_d               gmx_simd4_reduce_d_avx_128_fma
/* SIMD4 float/double conversion */
#define gmx_simd4_cvt_f2d                _mm256_cvtps_pd
#define gmx_simd4_cvt_d2f                _mm256_cvtpd_ps



/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_avx_128_fma(const float *        base,
                                              const gmx_int32_t    offset[],
                                              gmx_simd_float_t    &v0,
                                              gmx_simd_float_t    &v1,
                                              gmx_simd_float_t    &v2)
{
    __m128  t3;
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);

    assert((size_t)offset % 16 == 0);

    if ( (align & 0x3) == 0)
    {
        v0  = _mm_loadu_ps( base + align * offset[0] );
        v1  = _mm_loadu_ps( base + align * offset[1] );
        v2  = _mm_loadu_ps( base + align * offset[2] );
        t3  = _mm_loadu_ps( base + align * offset[3] );
    }
    else
    {
        v0  = gmx_mm_maskload_ps(base + align * offset[0], mask);
        v1  = gmx_mm_maskload_ps(base + align * offset[1], mask);
        v2  = gmx_mm_maskload_ps(base + align * offset[2], mask);
        t3  = gmx_mm_maskload_ps(base + align * offset[3], mask);
    }
    _MM_TRANSPOSE4_PS(v0, v1, v2, t3);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_avx_128_fma(float *              base,
                                                const gmx_int32_t    offset[],
                                                gmx_simd_float_t     v0,
                                                gmx_simd_float_t     v1,
                                                gmx_simd_float_t     v2)
{
    assert((size_t)offset % 16 == 0);

    __m128  t3   = _mm_setzero_ps();
    __m128i mask = _mm_set_epi32(0, -1, -1, -1);
    _MM_TRANSPOSE4_PS(v0, v1, v2, t3);
    gmx_mm_maskstore_ps(base + align * offset[0], mask, v0);
    gmx_mm_maskstore_ps(base + align * offset[1], mask, v1);
    gmx_mm_maskstore_ps(base + align * offset[2], mask, v2);
    gmx_mm_maskstore_ps(base + align * offset[3], mask, t3);

}

/* In the old group kernels we found it more efficient to transpose the data to store rather
 * than using maskload and maskstore. It might be worth to test again, but for now we assume
 * this is still the case, and rely on those version inherited from the SSE2 header.
 *
 * It is also worth testing if changing _mm_shuffle_ps() to _mm_permute_ps() could improve
 * throughput just-so-slightly.
 */

static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_f_avx_128_fma(gmx_simd_float_t   scalar,
                                                  gmx_simd_float_t  &triplets0,
                                                  gmx_simd_float_t  &triplets1,
                                                  gmx_simd_float_t  &triplets2)
{
    triplets0 = _mm_permute_ps(scalar, _MM_SHUFFLE(1, 0, 0, 0));
    triplets1 = _mm_permute_ps(scalar, _MM_SHUFFLE(2, 2, 1, 1));
    triplets2 = _mm_permute_ps(scalar, _MM_SHUFFLE(3, 3, 3, 2));
}

/* Permute might be a small help for gmx_simd_reduce_incr_4_return_sum_f(). Test in kernels later. */
#endif /* __cplusplus */




/****************************************************
 * Double precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd_expand_scalars_to_triplets_d_avx_128_fma(gmx_simd_double_t   scalar,
                                                  gmx_simd_double_t  &triplets0,
                                                  gmx_simd_double_t  &triplets1,
                                                  gmx_simd_double_t  &triplets2)
{
    triplets0 = _mm_permute_pd(scalar, _MM_SHUFFLE2(0, 0));
    triplets1 = _mm_permute_pd(scalar, _MM_SHUFFLE2(1, 0));
    triplets2 = _mm_permute_pd(scalar, _MM_SHUFFLE2(1, 1));
}


static gmx_inline double gmx_simdcall
gmx_simd4_reduce_d_avx_128_fma(__m256d a)
{
    double  f;
    __m128d a0, a1;
    /* test with shuffle & add as an alternative to hadd later */
    a  = _mm256_hadd_pd(a, a);
    a0 = _mm256_castpd256_pd128(a);
    a1 = _mm256_extractf128_pd(a, 0x1);
    a0 = _mm_add_sd(a0, a1);
    _mm_store_sd(&f, a0);
    return f;
}
#endif


static gmx_inline double gmx_simdcall
gmx_simd4_dotproduct3_d_avx_128_fma(__m256d a, __m256d b)
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

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_d_avx_128_fma(gmx_simd4_double_t &v0, gmx_simd4_double_t &v1,
                                  gmx_simd4_double_t &v2, gmx_simd4_double_t &v3)
{
    __m256d t1, t2, t3, t4;
    t1  = _mm256_unpacklo_pd(v0, v1);
    t2  = _mm256_unpackhi_pd(v0, v1);
    t3  = _mm256_unpacklo_pd(v2, v3);
    t4  = _mm256_unpackhi_pd(v2, v3);
    v0  = _mm256_permute2f128_pd(t1, t3, 0x20);
    v1  = _mm256_permute2f128_pd(t2, t4, 0x20);
    v2  = _mm256_permute2f128_pd(t1, t3, 0x31);
    v3  = _mm256_permute2f128_pd(t2, t4, 0x31);
}
#endif /* __cplusplus */

#endif /* GMX_SIMD_IMPL_X86_AVX_128_FMA_H */
