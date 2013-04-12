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
#ifndef _nbnxn_kernel_simd_utils_x86_256s_h_
#define _nbnxn_kernel_simd_utils_x86_256s_h_

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
static gmx_inline void
gmx_shuffle_4_ps_fil01_to_2_ps(__m128 in0, __m128 in1, __m128 in2, __m128 in3,
                               __m128 *out0, __m128 *out1)
{
    __m128 _c01, _c23;

    _c01  = _mm_movelh_ps(in0, in1);
    _c23  = _mm_movelh_ps(in2, in3);
    *out0 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));
    *out1 = _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(3, 1, 3, 1));
}

/* Collect element 2 of the 4 inputs to out */
static gmx_inline __m128
gmx_shuffle_4_ps_fil2_to_1_ps(__m128 in0, __m128 in1, __m128 in2, __m128 in3)
{
    __m128 _c01, _c23;

    _c01 = _mm_shuffle_ps(in0, in1, _MM_SHUFFLE(3, 2, 3, 2));
    _c23 = _mm_shuffle_ps(in2, in3, _MM_SHUFFLE(3, 2, 3, 2));

    return _mm_shuffle_ps(_c01, _c23, _MM_SHUFFLE(2, 0, 2, 0));
}

/* Sum the elements within each input register and store the sums in out */
static gmx_inline __m128
gmx_mm_transpose_sum4_pr(__m256 in0, __m256 in1,
                         __m256 in2, __m256 in3)
{
    in0 = _mm256_hadd_ps(in0, in1);
    in2 = _mm256_hadd_ps(in2, in3);
    in1 = _mm256_hadd_ps(in0, in2);

    return _mm_add_ps(_mm256_castps256_ps128(in1),
                      _mm256_extractf128_ps(in1, 1));
}

/* Sum the elements of halfs of each input register and store sums in out */
static gmx_inline __m128
gmx_mm_transpose_sum4h_pr(__m256 in0, __m256 in2)
{
    in0 = _mm256_hadd_ps(in0, _mm256_setzero_ps());
    in2 = _mm256_hadd_ps(in2, _mm256_setzero_ps());
    in0 = _mm256_hadd_ps(in0, in2);
    in2 = _mm256_permute_ps(in0, _MM_SHUFFLE(2, 3, 0, 1));

    return _mm_add_ps(_mm256_castps256_ps128(in0), _mm256_extractf128_ps(in2, 1));
}

/* Put two 128-bit 4-float registers into one 256-bit 8-float register */
static gmx_inline __m256
gmx_2_mm_to_m256(__m128 in0, __m128 in1)
{
    return _mm256_insertf128_ps(_mm256_castps128_ps256(in0), in1, 1);
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    __m256 *c6_S, __m256 *c12_S)
{
    __m128 clj_S[UNROLLJ], c6t_S[2], c12t_S[2];
    int    p;

    for (p = 0; p < UNROLLJ; p++)
    {
        /* Here we load 4 aligned floats, but we need just 2 */
        clj_S[p] = _mm_load_ps(nbfp+type[aj+p]*NBFP_STRIDE);
    }
    gmx_shuffle_4_ps_fil01_to_2_ps(clj_S[0], clj_S[1], clj_S[2], clj_S[3],
                                   &c6t_S[0], &c12t_S[0]);
    gmx_shuffle_4_ps_fil01_to_2_ps(clj_S[4], clj_S[5], clj_S[6], clj_S[7],
                                   &c6t_S[1], &c12t_S[1]);

    *c6_S  = gmx_2_mm_to_m256(c6t_S[0], c6t_S[1]);
    *c12_S = gmx_2_mm_to_m256(c12t_S[0], c12t_S[1]);
}

static gmx_inline void
load_lj_pair_params2(const real *nbfp0, const real *nbfp1,
                     const int *type, int aj,
                     __m256 *c6_S, __m256 *c12_S)
{
    __m128 clj_S0[UNROLLJ], clj_S1[UNROLLJ], c6t_S[2], c12t_S[2];
    int    p;

    for (p = 0; p < UNROLLJ; p++)
    {
        /* Here we load 4 aligned floats, but we need just 2 */
        clj_S0[p] = _mm_load_ps(nbfp0+type[aj+p]*NBFP_STRIDE);
    }
    for (p = 0; p < UNROLLJ; p++)
    {
        /* Here we load 4 aligned floats, but we need just 2 */
        clj_S1[p] = _mm_load_ps(nbfp1+type[aj+p]*NBFP_STRIDE);
    }
    gmx_shuffle_4_ps_fil01_to_2_ps(clj_S0[0], clj_S0[1], clj_S0[2], clj_S0[3],
                                   &c6t_S[0], &c12t_S[0]);
    gmx_shuffle_4_ps_fil01_to_2_ps(clj_S1[0], clj_S1[1], clj_S1[2], clj_S1[3],
                                   &c6t_S[1], &c12t_S[1]);

    *c6_S  = gmx_2_mm_to_m256(c6t_S[0], c6t_S[1]);
    *c12_S = gmx_2_mm_to_m256(c12t_S[0], c12t_S[1]);
}


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

static gmx_inline void
load_table_f(const real *tab_coul_FDV0, gmx_epi32 ti_S, int *ti,
             __m256 *ctab0_S, __m256 *ctab1_S)
{
    __m128 ctab_S[8], ctabt_S[4];
    int    j;

    /* Bit shifting would be faster, but AVX doesn't support that */
    _mm256_store_si256((__m256i *)ti, ti_S);
    for (j = 0; j < 8; j++)
    {
        ctab_S[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);
    }
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[0], ctab_S[1], ctab_S[2], ctab_S[3],
                                   &ctabt_S[0], &ctabt_S[2]);
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[4], ctab_S[5], ctab_S[6], ctab_S[7],
                                   &ctabt_S[1], &ctabt_S[3]);

    *ctab0_S = gmx_2_mm_to_m256(ctabt_S[0], ctabt_S[1]);
    *ctab1_S = gmx_2_mm_to_m256(ctabt_S[2], ctabt_S[3]);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_FDV0, gmx_epi32 ti_S, int *ti,
               __m256 *ctab0_S, __m256 *ctab1_S, __m256 *ctabv_S)
{
    __m128 ctab_S[8], ctabt_S[4], ctabvt_S[2];
    int    j;

    /* Bit shifting would be faster, but AVX doesn't support that */
    _mm256_store_si256((__m256i *)ti, ti_S);
    for (j = 0; j < 8; j++)
    {
        ctab_S[j] = _mm_load_ps(tab_coul_FDV0+ti[j]*4);
    }
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[0], ctab_S[1], ctab_S[2], ctab_S[3],
                                   &ctabt_S[0], &ctabt_S[2]);
    gmx_shuffle_4_ps_fil01_to_2_ps(ctab_S[4], ctab_S[5], ctab_S[6], ctab_S[7],
                                   &ctabt_S[1], &ctabt_S[3]);

    *ctab0_S = gmx_2_mm_to_m256(ctabt_S[0], ctabt_S[1]);
    *ctab1_S = gmx_2_mm_to_m256(ctabt_S[2], ctabt_S[3]);

    ctabvt_S[0] = gmx_shuffle_4_ps_fil2_to_1_ps(ctab_S[0], ctab_S[1],
                                                ctab_S[2], ctab_S[3]);
    ctabvt_S[1] = gmx_shuffle_4_ps_fil2_to_1_ps(ctab_S[4], ctab_S[5],
                                                ctab_S[6], ctab_S[7]);

    *ctabv_S = gmx_2_mm_to_m256(ctabvt_S[0], ctabvt_S[1]);
}


#ifdef gmx_mm_hpr
/* As add_ener_grp, but for two groups of UNROLLJ/2 stored in
 * a single SIMD register.
 */
static inline void
add_ener_grp_halves(__m256 e_S, real *v0, real *v1, const int *offset_jj)
{
    gmx_mm_hpr e_S0, e_S1;
    int        jj;

    e_S0 = _mm256_extractf128_ps(e_S, 0);
    e_S1 = _mm256_extractf128_ps(e_S, 1);

    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(v_S, v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_S, e_S0));
    }
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(v_S, v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_S, e_S1));
    }
}
#endif

#endif /* _nbnxn_kernel_simd_utils_x86_s256s_h_ */
