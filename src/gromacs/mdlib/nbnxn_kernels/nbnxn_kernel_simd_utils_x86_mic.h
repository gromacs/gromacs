/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifndef _nbnxn_kernel_simd_utils_x86_mic_h_
#define _nbnxn_kernel_simd_utils_x86_mic_h_

typedef gmx_simd_int32_t      gmx_exclfilter;
static const int filter_stride = GMX_SIMD_INT32_WIDTH/GMX_SIMD_REAL_WIDTH;

#define PERM_LOW2HIGH _MM_PERM_BABA
#define PERM_HIGH2LOW _MM_PERM_DCDC

#define mask_loh _mm512_int2mask(0x00FF) /* would be better a constant - but can't initialize with a function call. */
#define mask_hih _mm512_int2mask(0xFF00)

/* Half-width SIMD real type */
#ifdef GMX_SIMD_X86_MIC
typedef __m512 gmx_mm_hpr; /* high half is ignored */
typedef __m512i gmx_mm_hepi;
#else
typedef __m256 gmx_mm_hpr;
typedef __m256i gmx_mm_hepi;
#endif


/* Half-width SIMD operations */

/* Load reals at half-width aligned pointer b into half-width SIMD register a */
static gmx_inline void
gmx_load_hpr(gmx_mm_hpr *a, const real *b)
{
#ifdef GMX_SIMD_X86_MIC
    *a = _mm512_loadunpacklo_ps(_mm512_undefined_ps(), b);
#else
    *a = _mm256_load_ps(b);
#endif
}

/* Set all entries in half-width SIMD register *a to b */
static gmx_inline void
gmx_set1_hpr(gmx_mm_hpr *a, real b)
{
#ifdef GMX_SIMD_X86_MIC
    *a = _mm512_set1_ps(b);
#else
    *a = _mm256_set1_ps(b);
#endif
}

static gmx_inline void
gmx_2hpr_to_pr(gmx_mm_hpr a, gmx_mm_hpr b, gmx_simd_float_t *c);

/* Load one real at b and one real at b+1 into halves of a, respectively */
static gmx_inline void
gmx_load1p1_pr(gmx_simd_float_t *a, const real *b)
{
#ifdef GMX_SIMD_X86_MIC
    *a = _mm512_mask_extload_ps(_mm512_extload_ps(b, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE), mask_hih,
                                b+1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
#else
    gmx_2hpr_to_pr(_mm256_broadcast_ss(b), _mm256_broadcast_ss(b+1), a);
#endif
}

/* Load reals at half-width aligned pointer b into two halves of a */
static gmx_inline void
gmx_loaddh_pr(gmx_simd_float_t *a, const real *b)
{
#ifdef GMX_SIMD_X86_MIC
    *a = _mm512_permute4f128_ps(_mm512_loadunpacklo_ps(_mm512_undefined_ps(), b), PERM_LOW2HIGH);
#else
    __m256 t = _mm256_load_ps(b);
    gmx_2hpr_to_pr(t, t, a);    
#endif
}

/* Store half-width SIMD register b into half width aligned memory a */
static gmx_inline void
gmx_store_hpr(real *a, gmx_mm_hpr b)
{
#ifdef GMX_SIMD_X86_MIC
    _mm512_mask_packstorelo_ps(a, mask_loh, b);
#else
    _mm256_store_ps(a,b);
#endif
}

#ifdef GMX_SIMD_X86_MIC
#define gmx_add_hpr _mm512_add_ps
#define gmx_sub_hpr _mm512_sub_ps
#else
#define gmx_add_hpr _mm256_add_ps
#define gmx_sub_hpr _mm256_sub_ps
#endif

/* Sum over 4 half SIMD registers */
static gmx_inline gmx_mm_hpr
gmx_sum4_hpr(gmx_simd_float_t a, gmx_simd_float_t b)
{
    a = _mm512_add_ps(a, b);
    b = _mm512_permute4f128_ps(a, PERM_HIGH2LOW);
    a = _mm512_add_ps(a, b);
#ifdef GMX_SIMD_X86_MIC
    return a;
#else
    return _mm512_castps512_ps256(a);
#endif
}

#ifdef GMX_SIMD_X86_MIC
#define gmx_simd4_setr_f _mm512_setr4_ps
#else
#define gmx_simd4_setr_f _mm_setr_ps 
#endif

/* Sum the elements of halfs of each input register and store sums in out */
static gmx_inline gmx_simd4_real_t
gmx_mm_transpose_sum4h_pr(gmx_simd_float_t a, gmx_simd_float_t b)
{
    return gmx_simd4_setr_f(_mm512_mask_reduce_add_ps(mask_loh, a),
                            _mm512_mask_reduce_add_ps(mask_hih, a),
                            _mm512_mask_reduce_add_ps(mask_loh, b),
                            _mm512_mask_reduce_add_ps(mask_hih, b));
}

static gmx_inline void
gmx_pr_to_2hpr(gmx_simd_float_t a, gmx_mm_hpr *b, gmx_mm_hpr *c)
{
#ifdef GMX_SIMD_X86_MIC
    *b = a;
    *c = _mm512_permute4f128_ps(a, PERM_HIGH2LOW);
#else
/* *b = _mm512_extractf32x8_ps(a, 0); //TODO: does that actually produce an instruction? Shouldn't be required. Did it for consistency with avx256. But cast should be sufficient
   *c = _mm512_extractf32x8_ps(a, 1); //this compiles but the intrinsics guide says it is only DQ. What's correct? */
    *b = _mm512_castps512_ps256(a);
    *c = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(a), 1));
#endif
}

static gmx_inline void
gmx_2hpr_to_pr(gmx_mm_hpr a, gmx_mm_hpr b, gmx_simd_float_t *c)
{
#ifdef GMX_SIMD_X86_MIC
    *c = _mm512_mask_permute4f128_ps(a, mask_hih, b, PERM_LOW2HIGH);
#else
    *c = _mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_castps256_ps512(a)), _mm256_castps_pd(b), 0x1));  //the more proper 32x8 instruction is only available in DQ. Alternative one can use the MIC version
#endif
}

static gmx_inline void
gmx_2hepi_to_epi(gmx_mm_hepi a, gmx_mm_hepi b, gmx_simd_int32_t *c)
{
#ifdef GMX_SIMD_X86_MIC
    *c = _mm512_mask_permute4f128_epi32(a, mask_hih, b, PERM_LOW2HIGH);
#else
    *c = _mm512_inserti64x4(_mm512_castsi256_si512(a), b, 0x1);
#endif
}


static gmx_inline void
gmx_2hpr512_to_pr(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_float_t *c)
{
    *c = _mm512_mask_permute4f128_ps(a, mask_hih, b, PERM_LOW2HIGH);
}

/* recombine the 2 high half into c */
static gmx_inline void
gmx_2hpr512_high_to_pr(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_float_t *c)
{
    *c = _mm512_mask_permute4f128_ps(b, mask_loh, a, PERM_HIGH2LOW);
}

static gmx_inline void
gmx_2hepi512_to_epi(gmx_simd_int32_t a, gmx_simd_int32_t b, gmx_simd_int32_t *c)
{
    *c = _mm512_mask_permute4f128_epi32(a, mask_hih, b, PERM_LOW2HIGH);
}

/* recombine the 2 high half into c */
static gmx_inline void
gmx_2hepi512_high_to_epi(gmx_simd_int32_t a, gmx_simd_int32_t b, gmx_simd_int32_t *c)
{
    *c = _mm512_mask_permute4f128_epi32(b, mask_loh, a, PERM_HIGH2LOW);
}

/* Align a stack-based thread-local working array. work-array (currently) not used by load_table_f*/
static gmx_inline int *
prepare_table_load_buffer(const int *array)
{
    return NULL;
}

/* Using TAB_FDV0 is slower (for non-analytical PME). For the _mm512_i32gather_ps it helps
   to have the 16 elements in fewer cache lines. This is why it is faster to do F&D toghether
   and low/high half after each other, then simply doing a gather for tab_coul_F and tab_coul_F+1.
   The ording of the 16 elements doesn't matter, so it doesn't help to get FD sorted as odd/even
   instead of low/high.
 */
static gmx_inline void
load_table_f(const real *tab_coul_F, gmx_simd_int32_t ti_S, int *ti,
             gmx_simd_float_t *ctab0_S, gmx_simd_float_t *ctab1_S)
{
    __m512i idx;
    __m512i ti1 = _mm512_add_epi32(ti_S, _mm512_set1_epi32(1)); /* incr by 1 for tab1 */
    gmx_2hepi512_to_epi(ti_S, ti1, &idx);
    __m512  tmp1 = _mm512_i32gather_ps(idx, tab_coul_F, sizeof(float));
    gmx_2hepi512_high_to_epi(ti_S, ti1, &idx);
    __m512  tmp2 = _mm512_i32gather_ps(idx, tab_coul_F, sizeof(float));

    gmx_2hpr512_to_pr(tmp1, tmp2, ctab0_S);
    gmx_2hpr512_high_to_pr(tmp1, tmp2, ctab1_S);

    *ctab1_S  = gmx_simd_sub_r(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               gmx_simd_int32_t ti_S, int *ti,
               gmx_simd_float_t *ctab0_S, gmx_simd_float_t *ctab1_S,
               gmx_simd_float_t *ctabv_S)
{
    load_table_f(tab_coul_F, ti_S, ti, ctab0_S, ctab1_S);
    *ctabv_S = _mm512_i32gather_ps(ti_S, tab_coul_V, sizeof(float));
}

static gmx_inline __m512
gmx_mm_transpose_sum4_pr(gmx_simd_float_t in0, gmx_simd_float_t in1,
                         gmx_simd_float_t in2, gmx_simd_float_t in3)
{
    return _mm512_setr4_ps(_mm512_reduce_add_ps(in0),
                           _mm512_reduce_add_ps(in1),
                           _mm512_reduce_add_ps(in2),
                           _mm512_reduce_add_ps(in3));
}

static gmx_inline void
load_lj_pair_params2(const real *nbfp0, const real *nbfp1,
                     const int *type, int aj,
                     gmx_simd_float_t *c6_S, gmx_simd_float_t *c12_S)
{
    __m512i idx;
#ifdef GMX_SIMD_X86_MIC
    __m512i idx0, idx1;

    /* load all 8 unaligned requires 2 load. */
    idx0 = _mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), type+aj);
    idx0 = _mm512_loadunpackhi_epi32(idx0, type+aj+16);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(nbfp_stride));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1)); /* incr by 1 for c12 */
#else
    __m256i idx0, idx1;

    /* load all 8 unaligned requires 2 load. */
    idx0 = _mm256_loadu_si256((__m256i const*)(type+aj));

    idx0 = _mm256_mullo_epi32(idx0, _mm256_set1_epi32(nbfp_stride));
    idx1 = _mm256_add_epi32(idx0, _mm256_set1_epi32(1)); /* incr by 1 for c12 */
#endif

    gmx_2hepi_to_epi(idx0, idx1, &idx);
    __m512 tmp1 = _mm512_i32gather_ps(idx, nbfp0, sizeof(float));
    __m512 tmp2 = _mm512_i32gather_ps(idx, nbfp1, sizeof(float));

    gmx_2hpr512_to_pr(tmp1, tmp2, c6_S);
    gmx_2hpr512_high_to_pr(tmp1, tmp2, c12_S);
}

/* Code for handling loading exclusions and converting them into
   interactions. */
#define gmx_load1_exclfilter _mm512_set1_epi32
#define gmx_load_exclusion_filter _mm512_load_epi32
#define gmx_checkbitmask_pb _mm512_test_epi32_mask

#endif /* _nbnxn_kernel_simd_utils_ref_h_ */
