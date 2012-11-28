/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
 * Copyright (c) 2012, by the GROMACS development team, led by
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
#ifndef _nbnxn_kernel_sse_utils_h_
#define _nbnxn_kernel_sse_utils_h_

/* This files contains all functions/macros for the SSE/AVX kernels
 * which have explicit dependencies on the j-size / SIMD-width, which
 * can be 2 (SSE-double), 4 (SSE-single,AVX-double) or 8 (AVX-single).
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

#define GMX_MM_TRANSPOSE2_OP_PD(in0,in1,out0,out1)                      \
{                                                                       \
    out0 = _mm_shuffle_pd(in0,in1,_MM_SHUFFLE2(0,0));                   \
    out1 = _mm_shuffle_pd(in0,in1,_MM_SHUFFLE2(1,1));                   \
}

#if defined GMX_MM128_HERE || !defined GMX_DOUBLE
#define GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(in0,in1,in2,in3,out0,out1)    \
{                                                                       \
    __m128 _c01,_c23;                                                   \
    _c01 = _mm_shuffle_ps(in0,in1,_MM_SHUFFLE(1,0,1,0));                \
    _c23 = _mm_shuffle_ps(in2,in3,_MM_SHUFFLE(1,0,1,0));                \
    out0 = _mm_shuffle_ps(_c01,_c23,_MM_SHUFFLE(2,0,2,0));              \
    out1 = _mm_shuffle_ps(_c01,_c23,_MM_SHUFFLE(3,1,3,1));              \
}
#else
#define GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(in0,in1,in2,in3,out0,out1)    \
{                                                                       \
    __m256d _c01,_c23;                                                  \
    _c01 = _mm256_shuffle_pd(in0,in1,_MM_SHUFFLE(1,0,1,0));             \
    _c23 = _mm256_shuffle_pd(in2,in3,_MM_SHUFFLE(1,0,1,0));             \
    out0 = _mm256_shuffle_pd(_c01,_c23,_MM_SHUFFLE(2,0,2,0));           \
    out1 = _mm256_shuffle_pd(_c01,_c23,_MM_SHUFFLE(3,1,3,1));           \
}
#endif

#define GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(in0,in1,in2,in3,out)           \
{                                                                       \
    __m128 _c01,_c23;                                                   \
    _c01 = _mm_shuffle_ps(in0,in1,_MM_SHUFFLE(3,2,3,2));                \
    _c23 = _mm_shuffle_ps(in2,in3,_MM_SHUFFLE(3,2,3,2));                \
    out  = _mm_shuffle_ps(_c01,_c23,_MM_SHUFFLE(2,0,2,0));              \
}

#ifndef GMX_MM256_HERE
#ifndef GMX_DOUBLE
#define GMX_MM_TRANSPOSE_SUM4_PR(i_SSE0,i_SSE1,i_SSE2,i_SSE3,o_SSE)     \
{                                                                       \
    _MM_TRANSPOSE4_PS(i_SSE0,i_SSE1,i_SSE2,i_SSE3);                     \
    i_SSE0 = _mm_add_ps(i_SSE0,i_SSE1);                                 \
    i_SSE2 = _mm_add_ps(i_SSE2,i_SSE3);                                 \
    o_SSE  = _mm_add_ps(i_SSE0,i_SSE2);                                 \
}
#else
#define GMX_MM_TRANSPOSE_SUM2_PD(i_SSE0,i_SSE1,o_SSE)                   \
{                                                                       \
    GMX_MM_TRANSPOSE2_PD(i_SSE0,i_SSE1);                                \
    o_SSE  = _mm_add_pd(i_SSE0,i_SSE1);                                 \
}
#endif
#else
#ifndef GMX_DOUBLE
#define GMX_MM_TRANSPOSE_SUM4_PR(i_SSE0,i_SSE1,i_SSE2,i_SSE3,o_SSE)     \
{                                                                       \
    i_SSE0 = _mm256_hadd_ps(i_SSE0,i_SSE1);                             \
    i_SSE2 = _mm256_hadd_ps(i_SSE2,i_SSE3);                             \
    i_SSE1 = _mm256_hadd_ps(i_SSE0,i_SSE2);                             \
    o_SSE  = _mm_add_ps(_mm256_castps256_ps128(i_SSE1),_mm256_extractf128_ps(i_SSE1,1)); \
}
#else
#define GMX_MM_TRANSPOSE_SUM4_PR(i_SSE0,i_SSE1,i_SSE2,i_SSE3,o_SSE)     \
{                                                                       \
    i_SSE0 = _mm256_hadd_pd(i_SSE0,i_SSE1);                             \
    i_SSE2 = _mm256_hadd_pd(i_SSE2,i_SSE3);                             \
    o_SSE  = _mm256_add_pd(_mm256_permute2f128_pd(i_SSE0,i_SSE2,0x20),_mm256_permute2f128_pd(i_SSE0,i_SSE2,0x31)); \
}
#endif
#endif

#ifdef GMX_MM128_HERE

static inline __m128
gmx_mm128_invsqrt_ps_single(__m128 x)
{
    const __m128 half  = _mm_set_ps(0.5,0.5,0.5,0.5);
    const __m128 three = _mm_set_ps(3.0,3.0,3.0,3.0);
    
    __m128 lu = _mm_rsqrt_ps(x);
    
    return _mm_mul_ps(half,_mm_mul_ps(_mm_sub_ps(three,_mm_mul_ps(_mm_mul_ps(lu,lu),x)),lu));
}

/* Do 2/4 double precision invsqrt operations.
 * Doing the SSE rsqrt and the first Newton Raphson iteration
 * in single precision gives full double precision accuracy.
 * The speed is more than twice as fast as two gmx_mm_invsqrt_pd calls.
 */
#define GMX_MM128_INVSQRT2_PD(i_SSE0,i_SSE1,o_SSE0,o_SSE1)              \
{                                                                       \
    const __m128d half  = _mm_set1_pd(0.5);                             \
    const __m128d three = _mm_set1_pd(3.0);                             \
    __m128  s_SSE,ir_SSE;                                               \
    __m128d lu0,lu1;                                                    \
                                                                        \
    s_SSE  = _mm_movelh_ps(_mm_cvtpd_ps(i_SSE0),_mm_cvtpd_ps(i_SSE1));  \
    ir_SSE = gmx_mm128_invsqrt_ps_single(s_SSE);                        \
    lu0    = _mm_cvtps_pd(ir_SSE);                                      \
    lu1    = _mm_cvtps_pd(_mm_movehl_ps(ir_SSE,ir_SSE));                \
    o_SSE0 = _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu0,lu0),i_SSE0)),lu0)); \
    o_SSE1 = _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu1,lu1),i_SSE1)),lu1)); \
}

#define GMX_MM_INVSQRT2_PD GMX_MM128_INVSQRT2_PD

#endif

#ifdef GMX_MM256_HERE

static inline __m256
gmx_mm256_invsqrt_ps_single(__m256 x)
{
    const __m256 half  = _mm256_set_ps(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5);
    const __m256 three = _mm256_set_ps(3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0);
    
    __m256 lu = _mm256_rsqrt_ps(x);
    
    return _mm256_mul_ps(half,_mm256_mul_ps(_mm256_sub_ps(three,_mm256_mul_ps(_mm256_mul_ps(lu,lu),x)),lu));
}

#define GMX_MM256_INVSQRT2_PD(i_SSE0,i_SSE1,o_SSE0,o_SSE1)              \
{                                                                       \
    const __m256d half  = _mm256_set1_pd(0.5);                          \
    const __m256d three = _mm256_set1_pd(3.0);                          \
    __m256  s_SSE,ir_SSE;                                               \
    __m256d lu0,lu1;                                                    \
                                                                        \
    s_SSE  = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm256_cvtpd_ps(i_SSE0)),_mm256_cvtpd_ps(i_SSE1),1); \
    ir_SSE = gmx_mm256_invsqrt_ps_single(s_SSE);                        \
    lu0    = _mm256_cvtps_pd(_mm256_castps256_ps128(ir_SSE));           \
    lu1    = _mm256_cvtps_pd(_mm256_extractf128_ps(ir_SSE,1));          \
    o_SSE0 = _mm256_mul_pd(half,_mm256_mul_pd(_mm256_sub_pd(three,_mm256_mul_pd(_mm256_mul_pd(lu0,lu0),i_SSE0)),lu0)); \
    o_SSE1 = _mm256_mul_pd(half,_mm256_mul_pd(_mm256_sub_pd(three,_mm256_mul_pd(_mm256_mul_pd(lu1,lu1),i_SSE1)),lu1)); \
}

#define GMX_MM_INVSQRT2_PD GMX_MM256_INVSQRT2_PD

#endif

/* Force and energy table load and interpolation routines */

#if defined GMX_MM128_HERE && !defined GMX_DOUBLE

#define load_lj_pair_params(nbfp,type,aj,c6_SSE,c12_SSE)                \
{                                                                       \
    gmx_mm_pr clj_SSE[UNROLLJ];                                         \
    int p;                                                              \
                                                                        \
    for(p=0; p<UNROLLJ; p++)                                            \
    {                                                                   \
        /* Here we load 4 aligned floats, but we need just 2 */         \
        clj_SSE[p] = gmx_load_pr(nbfp+type[aj+p]*NBFP_STRIDE);          \
    }                                                                   \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE[0],clj_SSE[1],clj_SSE[2],clj_SSE[3],c6_SSE,c12_SSE); \
}

#endif

#if defined GMX_MM256_HERE && !defined GMX_DOUBLE

/* Put two 128-bit 4-float registers into one 256-bit 8-float register */
#define GMX_2_MM_TO_M256(in0,in1,out)                                   \
{                                                                       \
    out = _mm256_insertf128_ps(_mm256_castps128_ps256(in0),in1,1);      \
}

#define load_lj_pair_params(nbfp,type,aj,c6_SSE,c12_SSE)                \
{                                                                       \
    __m128 clj_SSE[UNROLLJ],c6t_SSE[2],c12t_SSE[2];                     \
    int p;                                                              \
                                                                        \
    for(p=0; p<UNROLLJ; p++)                                            \
    {                                                                   \
        /* Here we load 4 aligned floats, but we need just 2 */         \
        clj_SSE[p] = _mm_load_ps(nbfp+type[aj+p]*NBFP_STRIDE);          \
    }                                                                   \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE[0],clj_SSE[1],clj_SSE[2],clj_SSE[3],c6t_SSE[0],c12t_SSE[0]); \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(clj_SSE[4],clj_SSE[5],clj_SSE[6],clj_SSE[7],c6t_SSE[1],c12t_SSE[1]); \
                                                                        \
    GMX_2_MM_TO_M256(c6t_SSE[0],c6t_SSE[1],c6_SSE);                     \
    GMX_2_MM_TO_M256(c12t_SSE[0],c12t_SSE[1],c12_SSE);                  \
}

#endif

#if defined GMX_MM128_HERE && defined GMX_DOUBLE

#define load_lj_pair_params(nbfp,type,aj,c6_SSE,c12_SSE)                \
{                                                                       \
    gmx_mm_pr clj_SSE[UNROLLJ];                                         \
    int p;                                                              \
                                                                        \
    for(p=0; p<UNROLLJ; p++)                                            \
    {                                                                   \
        clj_SSE[p] = gmx_load_pr(nbfp+type[aj+p]*NBFP_STRIDE);          \
    }                                                                   \
    GMX_MM_TRANSPOSE2_OP_PD(clj_SSE[0],clj_SSE[1],c6_SSE,c12_SSE);      \
}

#endif

#if defined GMX_MM256_HERE && defined GMX_DOUBLE

#define load_lj_pair_params(nbfp,type,aj,c6_SSE,c12_SSE)                \
{                                                                       \
    __m128d clj_SSE[UNROLLJ],c6t_SSE[2],c12t_SSE[2];                    \
    int p;                                                              \
                                                                        \
    for(p=0; p<UNROLLJ; p++)                                            \
    {                                                                   \
        clj_SSE[p] = _mm_load_pd(nbfp+type[aj+p]*NBFP_STRIDE);          \
    }                                                                   \
    GMX_MM_TRANSPOSE2_OP_PD(clj_SSE[0],clj_SSE[1],c6t_SSE[0],c12t_SSE[0]); \
    GMX_MM_TRANSPOSE2_OP_PD(clj_SSE[2],clj_SSE[3],c6t_SSE[1],c12t_SSE[1]); \
    GMX_2_M128D_TO_M256D(c6t_SSE[0],c6t_SSE[1],c6_SSE);                 \
    GMX_2_M128D_TO_M256D(c12t_SSE[0],c12t_SSE[1],c12_SSE);              \
}

#endif


/* The load_table functions below are performance critical.
 * The routines issue UNROLLI*UNROLLJ _mm_load_ps calls.
 * As these all have latencies, scheduling is crucial.
 * The Intel compilers and CPUs seem to do a good job at this.
 * But AMD CPUs perform significantly worse with gcc than with icc.
 * Performance is improved a bit by using the extract function UNROLLJ times,
 * instead of doing an _mm_store_si128 for every i-particle.
 * With AVX this significantly deteriorates performance (8 extracts iso 4).
 * Because of this, the load_table_f macro always takes the ti parameter,
 * but it is only used with AVX.
 */

#if defined GMX_MM128_HERE && !defined GMX_DOUBLE

#define load_table_f(tab_coul_FDV0, ti_SSE, ti, ctab0_SSE, ctab1_SSE)   \
{                                                                       \
    int    idx[4];                                                      \
    __m128 ctab_SSE[4];                                                 \
                                                                        \
    /* Table has 4 entries, left-shift index by 2 */                    \
    ti_SSE = _mm_slli_epi32(ti_SSE,2);                                  \
    /* Without SSE4.1 the extract macro needs an immediate: unroll */   \
    idx[0] = gmx_mm_extract_epi32(ti_SSE,0);                            \
    ctab_SSE[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);                    \
    idx[1] = gmx_mm_extract_epi32(ti_SSE,1);                            \
    ctab_SSE[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);                    \
    idx[2] = gmx_mm_extract_epi32(ti_SSE,2);                            \
    ctab_SSE[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);                    \
    idx[3] = gmx_mm_extract_epi32(ti_SSE,3);                            \
    ctab_SSE[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);                    \
                                                                        \
    /* Shuffle the force table entries to a convenient order */         \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctab0_SSE,ctab1_SSE); \
}

#define load_table_f_v(tab_coul_FDV0, ti_SSE, ti, ctab0_SSE, ctab1_SSE, ctabv_SSE) \
{                                                                       \
    int    idx[4];                                                      \
    __m128 ctab_SSE[4];                                                 \
                                                                        \
    /* Table has 4 entries, left-shift index by 2 */                    \
    ti_SSE = _mm_slli_epi32(ti_SSE,2);                                  \
    /* Without SSE4.1 the extract macro needs an immediate: unroll */   \
    idx[0] = gmx_mm_extract_epi32(ti_SSE,0);                            \
    ctab_SSE[0] = _mm_load_ps(tab_coul_FDV0+idx[0]);                    \
    idx[1] = gmx_mm_extract_epi32(ti_SSE,1);                            \
    ctab_SSE[1] = _mm_load_ps(tab_coul_FDV0+idx[1]);                    \
    idx[2] = gmx_mm_extract_epi32(ti_SSE,2);                            \
    ctab_SSE[2] = _mm_load_ps(tab_coul_FDV0+idx[2]);                    \
    idx[3] = gmx_mm_extract_epi32(ti_SSE,3);                            \
    ctab_SSE[3] = _mm_load_ps(tab_coul_FDV0+idx[3]);                    \
                                                                        \
    /* Shuffle the force  table entries to a convenient order */        \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctab0_SSE,ctab1_SSE); \
    /* Shuffle the energy table entries to a convenient order */        \
    GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctabv_SSE); \
}

#endif

#if defined GMX_MM256_HERE && !defined GMX_DOUBLE

#define load_table_f(tab_coul_FDV0, ti_SSE, ti, ctab0_SSE, ctab1_SSE)   \
{                                                                       \
    __m128 ctab_SSE[8],ctabt_SSE[4];                                    \
    int    j;                                                           \
                                                                        \
    /* Bit shifting would be faster, but AVX doesn't support that */    \
    _mm256_store_si256((__m256i *)ti,ti_SSE);                           \
    for(j=0; j<8; j++)                                                  \
    {                                                                   \
        ctab_SSE[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);               \
    }                                                                   \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctabt_SSE[0],ctabt_SSE[2]); \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[4],ctab_SSE[5],ctab_SSE[6],ctab_SSE[7],ctabt_SSE[1],ctabt_SSE[3]); \
                                                                        \
    GMX_2_MM_TO_M256(ctabt_SSE[0],ctabt_SSE[1],ctab0_SSE);              \
    GMX_2_MM_TO_M256(ctabt_SSE[2],ctabt_SSE[3],ctab1_SSE);              \
}

#define load_table_f_v(tab_coul_FDV0, ti_SSE, ti, ctab0_SSE, ctab1_SSE, ctabv_SSE) \
{                                                                       \
    __m128 ctab_SSE[8],ctabt_SSE[4],ctabvt_SSE[2];                      \
    int    j;                                                           \
                                                                        \
    /* Bit shifting would be faster, but AVX doesn't support that */    \
    _mm256_store_si256((__m256i *)ti,ti_SSE);                           \
    for(j=0; j<8; j++)                                                  \
    {                                                                   \
        ctab_SSE[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);               \
    }                                                                   \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctabt_SSE[0],ctabt_SSE[2]); \
    GMX_MM_SHUFFLE_4_PS_FIL01_TO_2_PS(ctab_SSE[4],ctab_SSE[5],ctab_SSE[6],ctab_SSE[7],ctabt_SSE[1],ctabt_SSE[3]); \
                                                                        \
    GMX_2_MM_TO_M256(ctabt_SSE[0],ctabt_SSE[1],ctab0_SSE);              \
    GMX_2_MM_TO_M256(ctabt_SSE[2],ctabt_SSE[3],ctab1_SSE);              \
                                                                        \
    GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE[0],ctab_SSE[1],ctab_SSE[2],ctab_SSE[3],ctabvt_SSE[0]); \
    GMX_MM_SHUFFLE_4_PS_FIL2_TO_1_PS(ctab_SSE[4],ctab_SSE[5],ctab_SSE[6],ctab_SSE[7],ctabvt_SSE[1]); \
                                                                        \
    GMX_2_MM_TO_M256(ctabvt_SSE[0],ctabvt_SSE[1],ctabv_SSE);            \
}

#endif

#if defined GMX_MM128_HERE && defined GMX_DOUBLE

#define load_table_f(tab_coul_F, ti_SSE, ti, ctab0_SSE, ctab1_SSE)      \
{                                                                       \
    int     idx[2];                                                     \
    __m128d ctab_SSE[2];                                                \
                                                                        \
    /* Without SSE4.1 the extract macro needs an immediate: unroll */   \
    idx[0] = gmx_mm_extract_epi32(ti_SSE,0);                            \
    ctab_SSE[0] = _mm_loadu_pd(tab_coul_F+idx[0]);                      \
    idx[1] = gmx_mm_extract_epi32(ti_SSE,1);                            \
    ctab_SSE[1] = _mm_loadu_pd(tab_coul_F+idx[1]);                      \
                                                                        \
    /* Shuffle the force table entries to a convenient order */         \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[0],ctab_SSE[1],ctab0_SSE,ctab1_SSE); \
    /* The second force table entry should contain the difference */    \
    ctab1_SSE = _mm_sub_pd(ctab1_SSE,ctab0_SSE);                        \
}

#define load_table_f_v(tab_coul_F, tab_coul_V, ti_SSE, ti, ctab0_SSE, ctab1_SSE, ctabv_SSE) \
{                                                                       \
    int     idx[2];                                                     \
    __m128d ctab_SSE[4];                                                \
                                                                        \
    /* Without SSE4.1 the extract macro needs an immediate: unroll */   \
    idx[0] = gmx_mm_extract_epi32(ti_SSE,0);                            \
    ctab_SSE[0] = _mm_loadu_pd(tab_coul_F+idx[0]);                      \
    idx[1] = gmx_mm_extract_epi32(ti_SSE,1);                            \
    ctab_SSE[1] = _mm_loadu_pd(tab_coul_F+idx[1]);                      \
                                                                        \
    /* Shuffle the force table entries to a convenient order */         \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[0],ctab_SSE[1],ctab0_SSE,ctab1_SSE); \
    /* The second force table entry should contain the difference */    \
    ctab1_SSE = _mm_sub_pd(ctab1_SSE,ctab0_SSE);                        \
                                                                        \
    ctab_SSE[2] = _mm_loadu_pd(tab_coul_V+idx[0]);                      \
    ctab_SSE[3] = _mm_loadu_pd(tab_coul_V+idx[1]);                      \
                                                                        \
    /* Shuffle the energy table entries to a single register */         \
    ctabv_SSE = _mm_shuffle_pd(ctab_SSE[2],ctab_SSE[3],_MM_SHUFFLE2(0,0)); \
}

#endif

#if defined GMX_MM256_HERE && defined GMX_DOUBLE

/* Put two 128-bit 2-double registers into one 256-bit 4-ouble register */
#define GMX_2_M128D_TO_M256D(in0,in1,out)                               \
{                                                                       \
    out = _mm256_insertf128_pd(_mm256_castpd128_pd256(in0),in1,1);      \
}

#define load_table_f(tab_coul_F, ti_SSE, ti, ctab0_SSE, ctab1_SSE)      \
{                                                                       \
    __m128d ctab_SSE[4],tr_SSE[4];                                      \
    int     j;                                                          \
                                                                        \
    _mm_store_si128((__m128i *)ti,ti_SSE);                              \
    for(j=0; j<4; j++)                                                  \
    {                                                                   \
        ctab_SSE[j] = _mm_loadu_pd(tab_coul_F+ti[j]);                   \
    }                                                                   \
    /* Shuffle the force table entries to a convenient order */         \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[0],ctab_SSE[1],tr_SSE[0],tr_SSE[1]); \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[2],ctab_SSE[3],tr_SSE[2],tr_SSE[3]); \
    GMX_2_M128D_TO_M256D(tr_SSE[0],tr_SSE[2],ctab0_SSE);                \
    GMX_2_M128D_TO_M256D(tr_SSE[1],tr_SSE[3],ctab1_SSE);                \
    /* The second force table entry should contain the difference */    \
    ctab1_SSE = _mm256_sub_pd(ctab1_SSE,ctab0_SSE);                     \
}

#define load_table_f_v(tab_coul_F, tab_coul_V, ti_SSE, ti, ctab0_SSE, ctab1_SSE, ctabv_SSE) \
{                                                                       \
    __m128d ctab_SSE[8],tr_SSE[4];                                      \
    int     j;                                                          \
                                                                        \
    _mm_store_si128((__m128i *)ti,ti_SSE);                              \
    for(j=0; j<4; j++)                                                  \
    {                                                                   \
        ctab_SSE[j] = _mm_loadu_pd(tab_coul_F+ti[j]);                   \
    }                                                                   \
    /* Shuffle the force table entries to a convenient order */         \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[0],ctab_SSE[1],tr_SSE[0],tr_SSE[1]); \
    GMX_MM_TRANSPOSE2_OP_PD(ctab_SSE[2],ctab_SSE[3],tr_SSE[2],tr_SSE[3]); \
    GMX_2_M128D_TO_M256D(tr_SSE[0],tr_SSE[2],ctab0_SSE);                \
    GMX_2_M128D_TO_M256D(tr_SSE[1],tr_SSE[3],ctab1_SSE);                \
    /* The second force table entry should contain the difference */    \
    ctab1_SSE = _mm256_sub_pd(ctab1_SSE,ctab0_SSE);                     \
                                                                        \
    for(j=0; j<4; j++)                                                  \
    {                                                                   \
        ctab_SSE[4+j] = _mm_loadu_pd(tab_coul_V+ti[j]);                 \
    }                                                                   \
    /* Shuffle the energy table entries to a single register */         \
    GMX_2_M128D_TO_M256D(_mm_shuffle_pd(ctab_SSE[4],ctab_SSE[5],_MM_SHUFFLE2(0,0)),_mm_shuffle_pd(ctab_SSE[6],ctab_SSE[7],_MM_SHUFFLE2(0,0)),ctabv_SSE); \
}

#endif


/* Add energy register to possibly multiple terms in the energy array.
 * This function is the same for SSE/AVX single/double.
 */
static inline void add_ener_grp(gmx_mm_pr e_SSE,real *v,int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for(jj=0; jj<(UNROLLJ/2); jj++)
    {
        gmx_mm_pr v_SSE;

        v_SSE = gmx_load_pr(v+offset_jj[jj]+jj*UNROLLJ);
        gmx_store_pr(v+offset_jj[jj]+jj*UNROLLJ,gmx_add_pr(v_SSE,e_SSE));
    }
}

#endif /* _nbnxn_kernel_sse_utils_h_ */
