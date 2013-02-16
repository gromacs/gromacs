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
 * modify it under the terms of the GNU LeSIMDr General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * LeSIMDr General Public License for more details.
 *
 * You should have received a copy of the GNU LeSIMDr General Public
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
#ifndef _nbnxn_kernel_sse_utils_h_
#define _nbnxn_kernel_sse_utils_h_

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

#ifdef GMX_X86_SSE2

/* Transpose 2 double precision registers */
#define GMX_MM_TRANSPOSE2_OP_PD(in0, in1, out0, out1)                      \
    {                                                                       \
        out0 = _mm_unpacklo_pd(in0, in1);                                    \
        out1 = _mm_unpackhi_pd(in0, in1);                                    \
    }

#if defined GMX_MM128_HERE || !defined GMX_DOUBLE
/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
#define GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(in0, in1, in2, in3, out0, out1)    \
    {                                                                       \
        __m128 _c01, _c23;                                                   \
        _c01 = _mm_movelh_ps(in0, in1);                                      \
        _c23 = _mm_movelh_ps(in2, in3);                                      \
        out0 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));              \
        out1 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(3, 1, 3, 1));              \
    }
#else
/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
#define GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(in0, in1, in2, in3, out0, out1)    \
    {                                                                       \
        __m256d _c01, _c23;                                                  \
        _c01 = _mm256_shuffle_pd(in0, in1, _MM_SHUFFLE(1, 0, 1, 0));             \
        _c23 = _mm256_shuffle_pd(in2, in3, _MM_SHUFFLE(1, 0, 1, 0));             \
        out0 = _mm256_shuffle_pd(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));           \
        out1 = _mm256_shuffle_pd(_c01, _c23, _MM_SHUFFLE(3, 1, 3, 1));           \
    }
#endif

/* Collect element 2 of the 4 inputs to out */
#define GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(in0, in1, in2, in3, out)           \
    {                                                                       \
        __m128 _c01, _c23;                                                   \
        _c01 = _mm_shuffle_ps(in0, in1, _MM_SHUFFLE(3, 2, 3, 2));                \
        _c23 = _mm_shuffle_ps(in2, in3, _MM_SHUFFLE(3, 2, 3, 2));                \
        out  = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));              \
    }

#ifndef GMX_MM256_HERE
#ifndef GMX_DOUBLE
/* Sum the elements within each input register and store the sums in out */
#define GMX_MM_TRANSPOSE_SUM4_PR(in0, in1, in2, in3, out)                   \
    {                                                                       \
        _MM_TRANSPOSE4_PS(in0, in1, in2, in3);                                 \
        in0  = _mm_add_ps(in0, in1);                                          \
        in2  = _mm_add_ps(in2, in3);                                          \
        out  = _mm_add_ps(in0, in2);                                         \
    }
#else
/* Sum the elements within each input register and store the sums in out */
#define GMX_MM_TRANSPOSE_SUM2_PD(in0, in1, out)                           \
    {                                                                       \
        GMX_MM_TRANSPOSE2_PD(in0, in1);                                      \
        out  = _mm_add_pd(in0, in1);                                         \
    }
#endif
#else
#ifndef GMX_DOUBLE
/* Sum the elements within each input register and store the sums in out */
#define GMX_MM_TRANSPOSE_SUM4_PR(in0, in1, in2, in3, out)                   \
    {                                                                       \
        in0 = _mm256_hadd_ps(in0, in1);                                      \
        in2 = _mm256_hadd_ps(in2, in3);                                      \
        in1 = _mm256_hadd_ps(in0, in2);                                      \
        out = _mm_add_ps(_mm256_castps256_ps128(in1), _mm256_extractf128_ps(in1, 1)); \
    }
/* Sum the elements of halfs of each input register and store sums in out */
#define GMX_MM_TRANSPOSE_SUM4H_PR(in0, in2, out)                          \
    {                                                                       \
        in0 = _mm256_hadd_ps(in0, _mm256_setzero_ps());                      \
        in2 = _mm256_hadd_ps(in2, _mm256_setzero_ps());                      \
        in0 = _mm256_hadd_ps(in0, in2);                                      \
        in2 = _mm256_permute_ps(in0, _MM_SHUFFLE(2, 3, 0, 1));                  \
        out = _mm_add_ps(_mm256_castps256_ps128(in0), _mm256_extractf128_ps(in2, 1)); \
    }
#else
/* Sum the elements within each input register and store the sums in out */
#define GMX_MM_TRANSPOSE_SUM4_PR(in0, in1, in2, in3, out)                   \
    {                                                                       \
        in0 = _mm256_hadd_pd(in0, in1);                                      \
        in2 = _mm256_hadd_pd(in2, in3);                                      \
        out = _mm256_add_pd(_mm256_permute2f128_pd(in0, in2, 0x20), _mm256_permute2f128_pd(in0, in2, 0x31)); \
    }
#endif
#endif

#ifdef GMX_MM128_HERE

static inline __m128
gmx_mm128_invsqrt_ps_single(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5, 0.5, 0.5, 0.5);
    const __m128 three = _mm_set_ps(3.0, 3.0, 3.0, 3.0);

    __m128       lu = _mm_rsqrt_ps(x);

    return _mm_mul_ps(half, _mm_mul_ps(_mm_sub_ps(three, _mm_mul_ps(_mm_mul_ps(lu, lu), x)), lu));
}

/* Do 2 double precision invsqrt operations.
 * Doing the SIMD rsqrt and the first Newton Raphson iteration
 * in single precision gives full double precision accuracy.
 * The speed is more than double that of two gmx_mm_invsqrt_pd calls.
 */
#define GMX_MM128_INVSQRT2_PD(in0, in1, out0, out1)                        \
    {                                                                       \
        const __m128d half  = _mm_set1_pd(0.5);                             \
        const __m128d three = _mm_set1_pd(3.0);                             \
        __m128        s, ir;                                                       \
        __m128d       lu0, lu1;                                                    \
                                                                        \
        s    = _mm_movelh_ps(_mm_cvtpd_ps(in0), _mm_cvtpd_ps(in1));          \
        ir   = gmx_mm128_invsqrt_ps_single(s);                              \
        lu0  = _mm_cvtps_pd(ir);                                            \
        lu1  = _mm_cvtps_pd(_mm_movehl_ps(ir, ir));                          \
        out0 = _mm_mul_pd(half, _mm_mul_pd(_mm_sub_pd(three, _mm_mul_pd(_mm_mul_pd(lu0, lu0), in0)), lu0)); \
        out1 = _mm_mul_pd(half, _mm_mul_pd(_mm_sub_pd(three, _mm_mul_pd(_mm_mul_pd(lu1, lu1), in1)), lu1)); \
    }

#define GMX_MM_INVSQRT2_PD GMX_MM128_INVSQRT2_PD

#endif

#ifdef GMX_MM256_HERE

static inline __m256
gmx_mm256_invsqrt_ps_single(__m256 x)
{
    const __m256 half  = _mm256_set_ps(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
    const __m256 three = _mm256_set_ps(3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0);

    __m256       lu = _mm256_rsqrt_ps(x);

    return _mm256_mul_ps(half, _mm256_mul_ps(_mm256_sub_ps(three, _mm256_mul_ps(_mm256_mul_ps(lu, lu), x)), lu));
}

/* Do 4 double precision invsqrt operations.
 * Doing the SIMD rsqrt and the first Newton Raphson iteration
 * in single precision gives full double precision accuracy.
 */
#define GMX_MM256_INVSQRT2_PD(in0, in1, out0, out1)                        \
    {                                                                       \
        const __m256d half  = _mm256_set1_pd(0.5);                          \
        const __m256d three = _mm256_set1_pd(3.0);                          \
        __m256        s, ir;                                                       \
        __m256d       lu0, lu1;                                                    \
                                                                        \
        s    = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm256_cvtpd_ps(in0)), _mm256_cvtpd_ps(in1), 1); \
        ir   = gmx_mm256_invsqrt_ps_single(s);                              \
        lu0  = _mm256_cvtps_pd(_mm256_castps256_ps128(ir));                 \
        lu1  = _mm256_cvtps_pd(_mm256_extractf128_ps(ir, 1));                \
        out0 = _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu0, lu0), in0)), lu0)); \
        out1 = _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu1, lu1), in1)), lu1)); \
    }

#define GMX_MM_INVSQRT2_PD GMX_MM256_INVSQRT2_PD

#endif

/* Force and energy table load and interpolation routines */

#if defined GMX_MM128_HERE && !defined GMX_DOUBLE

#define load_lj_pair_params(nbfp, type, aj, c6_SIMD, c12_SIMD)                \
    {                                                                       \
        gmx_mm_pr clj_SIMD[UNROLLJ];                                         \
        int       p;                                                              \
                                                                        \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            /* Here we load 4 aligned floats, but we need just 2 */         \
            clj_SIMD[p] = gmx_load_pr(nbfp+type[aj+p]*NBFP_STRIDE);          \
        }                                                                   \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SIMD[0], clj_SIMD[1], clj_SIMD[2], clj_SIMD[3], c6_SIMD, c12_SIMD); \
    }

#endif

#if defined GMX_MM256_HERE && !defined GMX_DOUBLE

/* Put two 128-bit 4-float registers into one 256-bit 8-float register */
#define GMX_2_MM_TO_M256(in0, in1, out)                                   \
    {                                                                       \
        out = _mm256_insertf128_ps(_mm256_castps128_ps256(in0), in1, 1);      \
    }

#define load_lj_pair_params(nbfp, type, aj, c6_SIMD, c12_SIMD)                \
    {                                                                       \
        __m128 clj_SIMD[UNROLLJ], c6t_SIMD[2], c12t_SIMD[2];                     \
        int    p;                                                              \
                                                                        \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            /* Here we load 4 aligned floats, but we need just 2 */         \
            clj_SIMD[p] = _mm_load_ps(nbfp+type[aj+p]*NBFP_STRIDE);          \
        }                                                                   \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SIMD[0], clj_SIMD[1], clj_SIMD[2], clj_SIMD[3], c6t_SIMD[0], c12t_SIMD[0]); \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SIMD[4], clj_SIMD[5], clj_SIMD[6], clj_SIMD[7], c6t_SIMD[1], c12t_SIMD[1]); \
                                                                        \
        GMX_2_MM_TO_M256(c6t_SIMD[0], c6t_SIMD[1], c6_SIMD);                     \
        GMX_2_MM_TO_M256(c12t_SIMD[0], c12t_SIMD[1], c12_SIMD);                  \
    }

#define load_lj_pair_params2(nbfp0, nbfp1, type, aj, c6_SIMD, c12_SIMD)        \
    {                                                                       \
        __m128 clj_SIMD0[UNROLLJ], clj_SIMD1[UNROLLJ], c6t_SIMD[2], c12t_SIMD[2];  \
        int    p;                                                              \
                                                                        \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            /* Here we load 4 aligned floats, but we need just 2 */         \
            clj_SIMD0[p] = _mm_load_ps(nbfp0+type[aj+p]*NBFP_STRIDE);        \
        }                                                                   \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            /* Here we load 4 aligned floats, but we need just 2 */         \
            clj_SIMD1[p] = _mm_load_ps(nbfp1+type[aj+p]*NBFP_STRIDE);        \
        }                                                                   \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SIMD0[0], clj_SIMD0[1], clj_SIMD0[2], clj_SIMD0[3], c6t_SIMD[0], c12t_SIMD[0]); \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SIMD1[0], clj_SIMD1[1], clj_SIMD1[2], clj_SIMD1[3], c6t_SIMD[1], c12t_SIMD[1]); \
                                                                        \
        GMX_2_MM_TO_M256(c6t_SIMD[0], c6t_SIMD[1], c6_SIMD);                     \
        GMX_2_MM_TO_M256(c12t_SIMD[0], c12t_SIMD[1], c12_SIMD);                  \
    }

#endif

#if defined GMX_MM128_HERE && defined GMX_DOUBLE

#define load_lj_pair_params(nbfp, type, aj, c6_SIMD, c12_SIMD)                \
    {                                                                       \
        gmx_mm_pr clj_SIMD[UNROLLJ];                                         \
        int       p;                                                              \
                                                                        \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            clj_SIMD[p] = gmx_load_pr(nbfp+type[aj+p]*NBFP_STRIDE);          \
        }                                                                   \
        GMX_MM_TRANSPOSE2_OP_PD(clj_SIMD[0], clj_SIMD[1], c6_SIMD, c12_SIMD);      \
    }

#endif

#if defined GMX_MM256_HERE && defined GMX_DOUBLE

#define load_lj_pair_params(nbfp, type, aj, c6_SIMD, c12_SIMD)                \
    {                                                                       \
        __m128d clj_SIMD[UNROLLJ], c6t_SIMD[2], c12t_SIMD[2];                    \
        int     p;                                                              \
                                                                        \
        for (p = 0; p < UNROLLJ; p++)                                            \
        {                                                                   \
            clj_SIMD[p] = _mm_load_pd(nbfp+type[aj+p]*NBFP_STRIDE);          \
        }                                                                   \
        GMX_MM_TRANSPOSE2_OP_PD(clj_SIMD[0], clj_SIMD[1], c6t_SIMD[0], c12t_SIMD[0]); \
        GMX_MM_TRANSPOSE2_OP_PD(clj_SIMD[2], clj_SIMD[3], c6t_SIMD[1], c12t_SIMD[1]); \
        GMX_2_M128D_TO_M256D(c6t_SIMD[0], c6t_SIMD[1], c6_SIMD);                 \
        GMX_2_M128D_TO_M256D(c12t_SIMD[0], c12t_SIMD[1], c12_SIMD);              \
    }

#endif


/* The load_table functions below are performance critical.
 * The routines issue UNROLLI*UNROLLJ _mm_load_ps calls.
 * As these all have latencies, scheduling is crucial.
 * The Intel compilers and CPUs seem to do a good job at this.
 * But AMD CPUs perform significantly worse with gcc than with icc.
 * Performance is improved a bit by using the extract function UNROLLJ times,
 * instead of doing an _mm_store_si128 for every i-particle.
 * This is only faster when we use FDV0 formatted tables, where we also need
 * to multiple the index by 4, which can be done by a SIMD bit shift.
 * With single precision AVX, 8 extracts are much slower than 1 store.
 * Because of this, the load_table_f macro always takes the ti parameter,
 * but it is only used with AVX.
 */

#if defined GMX_MM128_HERE && !defined GMX_DOUBLE

#define load_table_f(tab_coul_FDV0, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD)   \
    {                                                                       \
        int    idx[4];                                                      \
        __m128 ctab_SIMD[4];                                                 \
                                                                        \
        /* Table has 4 entries, left-shift index by 2 */                    \
        ti_SIMD = _mm_slli_epi32(ti_SIMD, 2);                                  \
        /* Without SIMD4.1 the extract macro needs an immediate: unroll */   \
        idx[0]      = gmx_mm_extract_epi32(ti_SIMD, 0);                            \
        ctab_SIMD[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);                    \
        idx[1]      = gmx_mm_extract_epi32(ti_SIMD, 1);                            \
        ctab_SIMD[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);                    \
        idx[2]      = gmx_mm_extract_epi32(ti_SIMD, 2);                            \
        ctab_SIMD[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);                    \
        idx[3]      = gmx_mm_extract_epi32(ti_SIMD, 3);                            \
        ctab_SIMD[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);                    \
                                                                        \
        /* Shuffle the force table entries to a convenient order */         \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctab0_SIMD, ctab1_SIMD); \
    }

#define load_table_f_v(tab_coul_FDV0, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD, ctabv_SIMD) \
    {                                                                       \
        int    idx[4];                                                      \
        __m128 ctab_SIMD[4];                                                 \
                                                                        \
        /* Table has 4 entries, left-shift index by 2 */                    \
        ti_SIMD = _mm_slli_epi32(ti_SIMD, 2);                                  \
        /* Without SIMD4.1 the extract macro needs an immediate: unroll */   \
        idx[0]      = gmx_mm_extract_epi32(ti_SIMD, 0);                            \
        ctab_SIMD[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);                    \
        idx[1]      = gmx_mm_extract_epi32(ti_SIMD, 1);                            \
        ctab_SIMD[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);                    \
        idx[2]      = gmx_mm_extract_epi32(ti_SIMD, 2);                            \
        ctab_SIMD[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);                    \
        idx[3]      = gmx_mm_extract_epi32(ti_SIMD, 3);                            \
        ctab_SIMD[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);                    \
                                                                        \
        /* Shuffle the force  table entries to a convenient order */        \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctab0_SIMD, ctab1_SIMD); \
        /* Shuffle the energy table entries to a convenient order */        \
        GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctabv_SIMD); \
    }

#endif

#if defined GMX_MM256_HERE && !defined GMX_DOUBLE

#define load_table_f(tab_coul_FDV0, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD)   \
    {                                                                       \
        __m128 ctab_SIMD[8], ctabt_SIMD[4];                                    \
        int    j;                                                           \
                                                                        \
        /* Bit shifting would be faster, but AVX doesn't support that */    \
        _mm256_store_si256((__m256i *)ti, ti_SIMD);                           \
        for (j = 0; j < 8; j++)                                                  \
        {                                                                   \
            ctab_SIMD[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);               \
        }                                                                   \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctabt_SIMD[0], ctabt_SIMD[2]); \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[4], ctab_SIMD[5], ctab_SIMD[6], ctab_SIMD[7], ctabt_SIMD[1], ctabt_SIMD[3]); \
                                                                        \
        GMX_2_MM_TO_M256(ctabt_SIMD[0], ctabt_SIMD[1], ctab0_SIMD);              \
        GMX_2_MM_TO_M256(ctabt_SIMD[2], ctabt_SIMD[3], ctab1_SIMD);              \
    }

#define load_table_f_v(tab_coul_FDV0, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD, ctabv_SIMD) \
    {                                                                       \
        __m128 ctab_SIMD[8], ctabt_SIMD[4], ctabvt_SIMD[2];                      \
        int    j;                                                           \
                                                                        \
        /* Bit shifting would be faster, but AVX doesn't support that */    \
        _mm256_store_si256((__m256i *)ti, ti_SIMD);                           \
        for (j = 0; j < 8; j++)                                                  \
        {                                                                   \
            ctab_SIMD[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);               \
        }                                                                   \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctabt_SIMD[0], ctabt_SIMD[2]); \
        GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SIMD[4], ctab_SIMD[5], ctab_SIMD[6], ctab_SIMD[7], ctabt_SIMD[1], ctabt_SIMD[3]); \
                                                                        \
        GMX_2_MM_TO_M256(ctabt_SIMD[0], ctabt_SIMD[1], ctab0_SIMD);              \
        GMX_2_MM_TO_M256(ctabt_SIMD[2], ctabt_SIMD[3], ctab1_SIMD);              \
                                                                        \
        GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SIMD[0], ctab_SIMD[1], ctab_SIMD[2], ctab_SIMD[3], ctabvt_SIMD[0]); \
        GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SIMD[4], ctab_SIMD[5], ctab_SIMD[6], ctab_SIMD[7], ctabvt_SIMD[1]); \
                                                                        \
        GMX_2_MM_TO_M256(ctabvt_SIMD[0], ctabvt_SIMD[1], ctabv_SIMD);            \
    }

#endif

#if defined GMX_MM128_HERE && defined GMX_DOUBLE

#define load_table_f(tab_coul_F, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD)      \
    {                                                                       \
        int     idx[2];                                                     \
        __m128d ctab_SIMD[2];                                                \
                                                                        \
        /* Without SIMD4.1 the extract macro needs an immediate: unroll */   \
        idx[0]      = gmx_mm_extract_epi32(ti_SIMD, 0);                            \
        ctab_SIMD[0] = _mm_loadu_pd(tab_coul_F+idx[0]);                      \
        idx[1]      = gmx_mm_extract_epi32(ti_SIMD, 1);                            \
        ctab_SIMD[1] = _mm_loadu_pd(tab_coul_F+idx[1]);                      \
                                                                        \
        /* Shuffle the force table entries to a convenient order */         \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[0], ctab_SIMD[1], ctab0_SIMD, ctab1_SIMD); \
        /* The second force table entry should contain the difference */    \
        ctab1_SIMD = _mm_sub_pd(ctab1_SIMD, ctab0_SIMD);                        \
    }

#define load_table_f_v(tab_coul_F, tab_coul_V, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD, ctabv_SIMD) \
    {                                                                       \
        int     idx[2];                                                     \
        __m128d ctab_SIMD[4];                                                \
                                                                        \
        /* Without SIMD4.1 the extract macro needs an immediate: unroll */   \
        idx[0]      = gmx_mm_extract_epi32(ti_SIMD, 0);                            \
        ctab_SIMD[0] = _mm_loadu_pd(tab_coul_F+idx[0]);                      \
        idx[1]      = gmx_mm_extract_epi32(ti_SIMD, 1);                            \
        ctab_SIMD[1] = _mm_loadu_pd(tab_coul_F+idx[1]);                      \
                                                                        \
        /* Shuffle the force table entries to a convenient order */         \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[0], ctab_SIMD[1], ctab0_SIMD, ctab1_SIMD); \
        /* The second force table entry should contain the difference */    \
        ctab1_SIMD = _mm_sub_pd(ctab1_SIMD, ctab0_SIMD);                        \
                                                                        \
        ctab_SIMD[2] = _mm_loadu_pd(tab_coul_V+idx[0]);                      \
        ctab_SIMD[3] = _mm_loadu_pd(tab_coul_V+idx[1]);                      \
                                                                        \
        /* Shuffle the energy table entries to a single register */         \
        ctabv_SIMD = _mm_shuffle_pd(ctab_SIMD[2], ctab_SIMD[3], _MM_SHUFFLE2(0, 0)); \
    }

#endif

#if defined GMX_MM256_HERE && defined GMX_DOUBLE

/* Put two 128-bit 2-double registers into one 256-bit 4-ouble register */
#define GMX_2_M128D_TO_M256D(in0, in1, out)                               \
    {                                                                       \
        out = _mm256_insertf128_pd(_mm256_castpd128_pd256(in0), in1, 1);      \
    }

#define load_table_f(tab_coul_F, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD)      \
    {                                                                       \
        __m128d ctab_SIMD[4], tr_SIMD[4];                                      \
        int     j;                                                          \
                                                                        \
        _mm_store_si128((__m128i *)ti, ti_SIMD);                              \
        for (j = 0; j < 4; j++)                                                  \
        {                                                                   \
            ctab_SIMD[j] = _mm_loadu_pd(tab_coul_F+ti[j]);                   \
        }                                                                   \
        /* Shuffle the force table entries to a convenient order */         \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[0], ctab_SIMD[1], tr_SIMD[0], tr_SIMD[1]); \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[2], ctab_SIMD[3], tr_SIMD[2], tr_SIMD[3]); \
        GMX_2_M128D_TO_M256D(tr_SIMD[0], tr_SIMD[2], ctab0_SIMD);                \
        GMX_2_M128D_TO_M256D(tr_SIMD[1], tr_SIMD[3], ctab1_SIMD);                \
        /* The second force table entry should contain the difference */    \
        ctab1_SIMD = _mm256_sub_pd(ctab1_SIMD, ctab0_SIMD);                     \
    }

#define load_table_f_v(tab_coul_F, tab_coul_V, ti_SIMD, ti, ctab0_SIMD, ctab1_SIMD, ctabv_SIMD) \
    {                                                                       \
        __m128d ctab_SIMD[8], tr_SIMD[4];                                      \
        int     j;                                                          \
                                                                        \
        _mm_store_si128((__m128i *)ti, ti_SIMD);                              \
        for (j = 0; j < 4; j++)                                                  \
        {                                                                   \
            ctab_SIMD[j] = _mm_loadu_pd(tab_coul_F+ti[j]);                   \
        }                                                                   \
        /* Shuffle the force table entries to a convenient order */         \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[0], ctab_SIMD[1], tr_SIMD[0], tr_SIMD[1]); \
        GMX_MM_TRANSPOSE2_OP_PD(ctab_SIMD[2], ctab_SIMD[3], tr_SIMD[2], tr_SIMD[3]); \
        GMX_2_M128D_TO_M256D(tr_SIMD[0], tr_SIMD[2], ctab0_SIMD);                \
        GMX_2_M128D_TO_M256D(tr_SIMD[1], tr_SIMD[3], ctab1_SIMD);                \
        /* The second force table entry should contain the difference */    \
        ctab1_SIMD = _mm256_sub_pd(ctab1_SIMD, ctab0_SIMD);                     \
                                                                        \
        for (j = 0; j < 4; j++)                                                  \
        {                                                                   \
            ctab_SIMD[4+j] = _mm_loadu_pd(tab_coul_V+ti[j]);                 \
        }                                                                   \
        /* Shuffle the energy table entries to a single register */         \
        GMX_2_M128D_TO_M256D(_mm_shuffle_pd(ctab_SIMD[4], ctab_SIMD[5], _MM_SHUFFLE2(0, 0)), _mm_shuffle_pd(ctab_SIMD[6], ctab_SIMD[7], _MM_SHUFFLE2(0, 0)), ctabv_SIMD); \
    }

#endif


/* Add energy register to possibly multiple terms in the energy array.
 * This function is the same for SIMD/AVX single/double.
 */
static inline void add_ener_grp(gmx_mm_pr e_SIMD, real *v, const int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_pr v_SIMD;

        v_SIMD = gmx_load_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE);
        gmx_store_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE, gmx_add_pr(v_SIMD, e_SIMD));
    }
}

#if defined GMX_X86_AVX_256 && GMX_SIMD_WIDTH_HERE == 8
/* As add_ener_grp above, but for two groups of UNROLLJ/2 stored in
 * a single SIMD register.
 */
static inline void add_ener_grp_halves(gmx_mm_pr e_SIMD,
                                       real *v0, real *v1, const int *offset_jj)
{
    gmx_mm_hpr e_SIMD0, e_SIMD1;
    int        jj;

    e_SIMD0 = _mm256_extractf128_ps(e_SIMD, 0);
    e_SIMD1 = _mm256_extractf128_ps(e_SIMD, 1);

    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_SIMD;

        v_SIMD = gmx_load_hpr(v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_SIMD, e_SIMD0));
    }
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_SIMD;

        v_SIMD = gmx_load_hpr(v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_SIMD, e_SIMD1));
    }
}
#endif

#endif /* GMX_X86_SSE2 */

#endif /* _nbnxn_kernel_sse_utils_h_ */
