/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_H
#define GMX_SIMD_IMPL_X86_AVX2_256_H

#include "config.h"

#include <math.h>

#include <immintrin.h>

/* x86 256-bit AVX2 SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for details
 */

/* Inherit parts of AVX2_256 from AVX_256 */
#include "gromacs/simd/impl_x86_avx_256/impl_x86_avx_256.h"
/* Increment over AVX_256 capabilities */
#define GMX_SIMD_X86_AVX2_256_OR_HIGHER

/* Override some capability definitions for things added in AVX2 */
#define GMX_SIMD_HAVE_FMA
#define GMX_SIMD_HAVE_FINT32_LOGICAL     /* AVX2 adds 256-bit integer shifts */
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS /* AVX2 adds 256-bit integer +,-,*  */

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#undef  gmx_simd_fmadd_f
#define gmx_simd_fmadd_f           _mm256_fmadd_ps
#undef  gmx_simd_fmsub_f
#define gmx_simd_fmsub_f           _mm256_fmsub_ps
#undef  gmx_simd_fnmadd_f
#define gmx_simd_fnmadd_f          _mm256_fnmadd_ps
#undef  gmx_simd_fnmsub_f
#define gmx_simd_fnmsub_f          _mm256_fnmsub_ps
#undef  gmx_simd_get_exponent_f
#define gmx_simd_get_exponent_f    gmx_simd_get_exponent_f_avx2_256
#undef  gmx_simd_set_exponent_f
#define gmx_simd_set_exponent_f    gmx_simd_set_exponent_f_avx2_256
/* Previously undefined logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi           _mm256_slli_epi32
#define gmx_simd_srli_fi           _mm256_srli_epi32
#define gmx_simd_and_fi            _mm256_and_si256
#define gmx_simd_andnot_fi         _mm256_andnot_si256
#define gmx_simd_or_fi             _mm256_or_si256
#define gmx_simd_xor_fi            _mm256_xor_si256
/* Previously undefined arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi            _mm256_add_epi32
#define gmx_simd_sub_fi            _mm256_sub_epi32
#define gmx_simd_mul_fi            _mm256_mullo_epi32
/* Previously undefined boolean ops on gmx_simd_fint32_t */
#define gmx_simd_cmpeq_fi          _mm256_cmpeq_epi32
#define gmx_simd_cmplt_fi(a, b)     _mm256_cmpgt_epi32(b, a)
#define gmx_simd_and_fib           _mm256_and_si256
#define gmx_simd_or_fib            _mm256_or_si256
#define gmx_simd_anytrue_fib       _mm256_movemask_epi8
#define gmx_simd_blendzero_fi      _mm256_and_si256
#define gmx_simd_blendnotzero_fi(a, sel) _mm256_andnot_si256(sel, a)
#define gmx_simd_blendv_fi         _mm256_blendv_epi8


/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#undef  gmx_simd_fmadd_d
#define gmx_simd_fmadd_d           _mm256_fmadd_pd
#undef  gmx_simd_fmsub_d
#define gmx_simd_fmsub_d           _mm256_fmsub_pd
#undef  gmx_simd_fnmadd_d
#define gmx_simd_fnmadd_d          _mm256_fnmadd_pd
#undef  gmx_simd_fnmsub_d
#define gmx_simd_fnmsub_d          _mm256_fnmsub_pd
#undef  gmx_simd_get_exponent_d
#define gmx_simd_get_exponent_d    gmx_simd_get_exponent_d_avx2_256
#undef  gmx_simd_set_exponent_d
#define gmx_simd_set_exponent_d    gmx_simd_set_exponent_d_avx2_256
#undef  gmx_simd_cvt_db2dib
#define gmx_simd_cvt_db2dib        gmx_simd_cvt_db2dib_avx2_256
#undef  gmx_simd_cvt_dib2db
#define gmx_simd_cvt_dib2db        gmx_simd_cvt_dib2db_avx2_256

/****************************************************
 *      SIMD4 SINGLE PRECISION IMPLEMENTATION       *
 ****************************************************/
#undef  gmx_simd4_fmadd_f
#define gmx_simd4_fmadd_f          _mm_fmadd_ps
#undef  gmx_simd4_fmsub_f
#define gmx_simd4_fmsub_f          _mm_fmsub_ps
#undef  gmx_simd4_fnmadd_f
#define gmx_simd4_fnmadd_f         _mm_fnmadd_ps
#undef  gmx_simd4_fnmsub_f
#define gmx_simd4_fnmsub_f         _mm_fnmsub_ps

/* No need to update SIMD4 double, since those instructions
 * are aliased to the general SIMD double instructions above.
 */

/*********************************************************
 * SIMD SINGLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_get_exponent_f_avx2_256(gmx_simd_float_t x)
{
    const __m256  expmask      = _mm256_castsi256_ps(_mm256_set1_epi32(0x7F800000));
    const __m256i expbias      = _mm256_set1_epi32(127);
    __m256i       iexp;

    iexp = _mm256_castps_si256(_mm256_and_ps(x, expmask));
    iexp = _mm256_sub_epi32(_mm256_srli_epi32(iexp, 23), expbias);
    return _mm256_cvtepi32_ps(iexp);
}

static gmx_inline gmx_simd_float_t gmx_simdcall
gmx_simd_set_exponent_f_avx2_256(gmx_simd_float_t x)
{
    const __m256i  expbias      = _mm256_set1_epi32(127);
    __m256i        iexp         = _mm256_cvtps_epi32(x);

    iexp = _mm256_slli_epi32(_mm256_add_epi32(iexp, expbias), 23);
    return _mm256_castsi256_ps(iexp);
}

/*********************************************************
 * SIMD DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 *********************************************************/
static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_get_exponent_d_avx2_256(gmx_simd_double_t x)
{
    const __m256d  expmask      = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FF0000000000000LL));
    const __m256i  expbias      = _mm256_set1_epi64x(1023LL);
    __m256i        iexp;
    __m128i        iexp128;

    iexp = _mm256_castpd_si256(_mm256_and_pd(x, expmask));
    iexp = _mm256_sub_epi64(_mm256_srli_epi64(iexp, 52), expbias);
    iexp = _mm256_shuffle_epi32(iexp, _MM_SHUFFLE(3, 1, 2, 0));

    iexp128 = _mm256_extractf128_si256(iexp, 1);
    iexp128 = _mm_unpacklo_epi64(_mm256_castsi256_si128(iexp), iexp128);
    return _mm256_cvtepi32_pd(iexp128);
}

static gmx_inline gmx_simd_double_t gmx_simdcall
gmx_simd_set_exponent_d_avx2_256(gmx_simd_double_t x)
{
    const __m256i  expbias      = _mm256_set1_epi64x(1023LL);
    __m256i        iexp         = _mm256_cvtepi32_epi64(_mm256_cvtpd_epi32(x));

    iexp = _mm256_slli_epi64(_mm256_add_epi64(iexp, expbias), 52);
    return _mm256_castsi256_pd(iexp);
}

static gmx_inline gmx_simd_dibool_t gmx_simdcall
gmx_simd_cvt_db2dib_avx2_256(gmx_simd_dbool_t a)
{
    __m128i ia = _mm256_castsi256_si128(_mm256_castpd_si256(a));
    __m128i ib = _mm256_extractf128_si256(_mm256_castpd_si256(a), 0x1);

    ia = _mm_packs_epi32(ia, ib);

    return ia;
}

static gmx_inline gmx_simd_dbool_t gmx_simdcall
gmx_simd_cvt_dib2db_avx2_256(gmx_simd_dibool_t ia)
{
    __m128d lo = _mm_castsi128_pd(_mm_unpacklo_epi32(ia, ia));
    __m128d hi = _mm_castsi128_pd(_mm_unpackhi_epi32(ia, ia));

    return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 0x1);
}

#endif /* GMX_SIMD_IMPL_X86_AVX2_256_H */
