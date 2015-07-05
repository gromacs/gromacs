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

#ifndef GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx2_256_common.h"

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

#endif /* GMX_SIMD_IMPL_X86_AVX2_256_SIMD_FLOAT_H */
