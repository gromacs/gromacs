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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx_256_common.h"
#include "impl_x86_avx_256_simd_double.h"

/****************************************************
 *      DOUBLE PRECISION SIMD4 IMPLEMENTATION       *
 ****************************************************/
#define gmx_simd4_double_t          gmx_simd_double_t
#define gmx_simd4_load_d            gmx_simd_load_d
#define gmx_simd4_load1_d           gmx_simd_load1_d
#define gmx_simd4_set1_d            gmx_simd_set1_d
#define gmx_simd4_store_d           gmx_simd_store_d
#define gmx_simd4_loadu_d           gmx_simd_loadu_d
#define gmx_simd4_storeu_d          gmx_simd_storeu_d
#define gmx_simd4_setzero_d         gmx_simd_setzero_d
#define gmx_simd4_add_d             gmx_simd_add_d
#define gmx_simd4_sub_d             gmx_simd_sub_d
#define gmx_simd4_mul_d             gmx_simd_mul_d
#define gmx_simd4_fmadd_d           gmx_simd_fmadd_d
#define gmx_simd4_fmsub_d           gmx_simd_fmsub_d
#define gmx_simd4_fnmadd_d          gmx_simd_fnmadd_d
#define gmx_simd4_fnmsub_d          gmx_simd_fnmsub_d
#define gmx_simd4_and_d             gmx_simd_and_d
#define gmx_simd4_andnot_d          gmx_simd_andnot_d
#define gmx_simd4_or_d              gmx_simd_or_d
#define gmx_simd4_xor_d             gmx_simd_xor_d
#define gmx_simd4_rsqrt_d           gmx_simd_rsqrt_d
#define gmx_simd4_fabs_d            gmx_simd_fabs_d
#define gmx_simd4_fneg_d            gmx_simd_fneg_d
#define gmx_simd4_max_d             gmx_simd_max_d
#define gmx_simd4_min_d             gmx_simd_min_d
#define gmx_simd4_round_d           gmx_simd_round_d
#define gmx_simd4_trunc_d           gmx_simd_trunc_d
#define gmx_simd4_dotproduct3_d     gmx_simd4_dotproduct3_d_avx_256
#define gmx_simd4_dbool_t           gmx_simd_dbool_t
#define gmx_simd4_cmpeq_d           gmx_simd_cmpeq_d
#define gmx_simd4_cmplt_d           gmx_simd_cmplt_d
#define gmx_simd4_cmple_d           gmx_simd_cmple_d
#define gmx_simd4_and_db            gmx_simd_and_db
#define gmx_simd4_or_db             gmx_simd_or_db
#define gmx_simd4_anytrue_db        gmx_simd_anytrue_db
#define gmx_simd4_blendzero_d       gmx_simd_blendzero_d
#define gmx_simd4_blendnotzero_d    gmx_simd_blendnotzero_d
#define gmx_simd4_blendv_d          gmx_simd_blendv_d
#define gmx_simd4_reduce_d          gmx_simd_reduce_d
/* SIMD4 float/double conversion */
#define gmx_simd4_cvt_f2d           _mm256_cvtps_pd
#define gmx_simd4_cvt_d2f           _mm256_cvtpd_ps

/* Implementation helpers */
static gmx_inline double gmx_simdcall
gmx_simd4_dotproduct3_d_avx_256(__m256d a, __m256d b)
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


#endif /* GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H */
