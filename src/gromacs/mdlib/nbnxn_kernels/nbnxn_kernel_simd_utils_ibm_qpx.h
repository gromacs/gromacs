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
#ifndef _nbnxn_kernel_simd_utils_ibm_qpx_h_
#define _nbnxn_kernel_simd_utils_ibm_qpx_h_

typedef gmx_simd_real_t gmx_exclfilter;
static const int filter_stride = 1;

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

/* Collect all [0123] elements of the 4 inputs to out[0123], respectively */
static gmx_inline void
gmx_transpose_4_ps(gmx_simd_real_t a, gmx_simd_real_t b,
                   gmx_simd_real_t c, gmx_simd_real_t d,
                   gmx_simd_real_t *out0, gmx_simd_real_t *out1,
                   gmx_simd_real_t *out2, gmx_simd_real_t *out3)
{
    /* Prepare control vectors for swizzling. In its third input,
       vec_perm accepts indices into the effective 8-wide SIMD vector
       created by concatenating its first two inputs. Those indices
       map data from the input vectors to the output vector.

       vec_gpci() converts an octal literal of the indices into the
       correct form for vec_perm() to use. That form is an octal digit
       in bits 0-2 of the mantissa of each double. */
    gmx_simd_real_t p6420 = vec_gpci(06420);
    gmx_simd_real_t p7531 = vec_gpci(07531);

    /* Four-way swizzle (i.e. transpose) of vectors a = a0a1a2a3, etc. */
    gmx_simd_real_t b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_simd_real_t b3b1a3a1 = vec_perm(a, b, p7531);
    gmx_simd_real_t d2d0c2c0 = vec_perm(c, d, p6420);
    gmx_simd_real_t d3d1c3c1 = vec_perm(c, d, p7531);
    *out0 = vec_perm(d2d0c2c0, b2b0a2a0, p7531);
    *out1 = vec_perm(d3d1c3c1, b3b1a3a1, p7531);
    *out2 = vec_perm(d2d0c2c0, b2b0a2a0, p6420);
    *out3 = vec_perm(d3d1c3c1, b3b1a3a1, p6420);
}

/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
static gmx_inline void
gmx_shuffle_4_ps_fil01_to_2_ps(gmx_simd_real_t a, gmx_simd_real_t b,
                               gmx_simd_real_t c, gmx_simd_real_t d,
                               gmx_simd_real_t *out0, gmx_simd_real_t *out1)
{
    gmx_simd_real_t p6420 = vec_gpci(06420);
    gmx_simd_real_t p7531 = vec_gpci(07531);

    /* Partial four-way swizzle of vectors a = a0a1a2a3, etc. */
    gmx_simd_real_t b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_simd_real_t b3b1a3a1 = vec_perm(a, b, p7531);
    gmx_simd_real_t d2d0c2c0 = vec_perm(c, d, p6420);
    gmx_simd_real_t d3d1c3c1 = vec_perm(c, d, p7531);
    *out0 = vec_perm(d2d0c2c0, b2b0a2a0, p7531);
    *out1 = vec_perm(d3d1c3c1, b3b1a3a1, p7531);
}

/* Collect element 2 of the 4 inputs to out */
static gmx_inline gmx_simd_real_t
gmx_shuffle_4_ps_fil2_to_1_ps(gmx_simd_real_t a, gmx_simd_real_t b,
                              gmx_simd_real_t c, gmx_simd_real_t d)
{
    gmx_simd_real_t p6420 = vec_gpci(06420);

    /* Partial four-way swizzle of vectors a = a0a1a2a3, etc. */
    gmx_simd_real_t b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_simd_real_t d2d0c2c0 = vec_perm(c, d, p6420);
    return vec_perm(d2d0c2c0, b2b0a2a0, p6420);
}

#ifdef TAB_FDV0
/* Align a stack-based thread-local working array. Table loads on QPX
 * use the array, but most other implementations do not. */
static gmx_inline int *
prepare_table_load_buffer(int *array)
{
    return gmx_simd_align_i(array);
}

static gmx_inline void
load_table_f(const real *tab_coul_FDV0, gmx_simd_int32_t ti_S, int *ti,
             gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S)
{
#ifdef NDEBUG
    /* Just like 256-bit AVX, we need to use memory to get indices
       into integer registers efficiently. */
    vec_st(ti_S, 0, ti);
#else
    vec_sta(ti_S, 0, ti);
#endif

    /* Here we load 4 aligned reals, but we need just 2 elements of each */
    gmx_simd_real_t a = gmx_simd_load_r(tab_coul_FDV0 + ti[0] * nbfp_stride);
    gmx_simd_real_t b = gmx_simd_load_r(tab_coul_FDV0 + ti[1] * nbfp_stride);
    gmx_simd_real_t c = gmx_simd_load_r(tab_coul_FDV0 + ti[2] * nbfp_stride);
    gmx_simd_real_t d = gmx_simd_load_r(tab_coul_FDV0 + ti[3] * nbfp_stride);

    gmx_shuffle_4_ps_fil01_to_2_ps(a, b, c, d, ctab0_S, ctab1_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_FDV0,
               gmx_simd_int32_t ti_S, int *ti,
               gmx_simd_real_t *ctab0_S, gmx_simd_real_t *ctab1_S,
               gmx_simd_real_t *ctabv_S)
{
#ifdef NDEBUG
    /* Just like 256-bit AVX, we need to use memory to get indices
       into integer registers efficiently. */
    vec_st(ti_S, 0, ti);
#else
    vec_sta(ti_S, 0, ti);
#endif

    /* Here we load 4 aligned reals, but we need just 3 elements of each. */
    gmx_simd_real_t a = gmx_simd_load_r(tab_coul_FDV0 + ti[0] * nbfp_stride);
    gmx_simd_real_t b = gmx_simd_load_r(tab_coul_FDV0 + ti[1] * nbfp_stride);
    gmx_simd_real_t c = gmx_simd_load_r(tab_coul_FDV0 + ti[2] * nbfp_stride);
    gmx_simd_real_t d = gmx_simd_load_r(tab_coul_FDV0 + ti[3] * nbfp_stride);

    gmx_shuffle_4_ps_fil01_to_2_ps(a, b, c, d, ctab0_S, ctab1_S);
    *ctabv_S = gmx_shuffle_4_ps_fil2_to_1_ps(a, b, c, d);
}
#else

/* Not required for BlueGene/Q */

#endif

/* Sum the elements within each input register and store the sums in out.
 */
static gmx_inline gmx_simd_real_t
gmx_mm_transpose_sum4_pr(gmx_simd_real_t a, gmx_simd_real_t b,
                         gmx_simd_real_t c, gmx_simd_real_t d)
{
    gmx_simd_real_t a0b0c0d0, a1b1c1d1, a2b2c2d2, a3b3c3d3;
    gmx_transpose_4_ps(a, b, c, d,
                       &a0b0c0d0,
                       &a1b1c1d1,
                       &a2b2c2d2,
                       &a3b3c3d3);
    /* Now reduce the transposed vectors */
    gmx_simd_real_t sum01 = gmx_simd_add_r(a0b0c0d0, a1b1c1d1);
    gmx_simd_real_t sim23 = gmx_simd_add_r(a2b2c2d2, a3b3c3d3);
    return gmx_simd_add_r(sum01, sim23);
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    gmx_simd_real_t *c6_S, gmx_simd_real_t *c12_S)
{
    /* Here we load 4 aligned reals, but we need just 2 elemnts of each. */
    gmx_simd_real_t a = gmx_simd_load_r(nbfp + type[aj+0] * nbfp_stride);
    gmx_simd_real_t b = gmx_simd_load_r(nbfp + type[aj+1] * nbfp_stride);
    gmx_simd_real_t c = gmx_simd_load_r(nbfp + type[aj+2] * nbfp_stride);
    gmx_simd_real_t d = gmx_simd_load_r(nbfp + type[aj+3] * nbfp_stride);

    gmx_shuffle_4_ps_fil01_to_2_ps(a, b, c, d, c6_S, c12_S);
}

/* Define USE_FUNCTIONS_FOR_QPX to get the static inline functions
 * that seem to exhaust xlC 12.1 during kernel compilation */

#ifndef USE_FUNCTIONS_FOR_QPX

#define gmx_load_exclusion_filter(a) vec_ldia(0, (int *) a)
#define gmx_load_interaction_mask_pb(a, b) vec_ld(a, (real *) b)

#else /* USE_FUNCTIONS_FOR_QPX */

static gmx_inline gmx_exclfilter gmx_load_exclusion_filter(const unsigned *a)
{
#ifdef NDEBUG
    return vec_ldia(0, (int *) a);
#else
    return vec_ldiaa(0, (int *) a);
#endif
}

/* Code for handling loading and applying exclusion masks. Note that
   parameter a is not treated like an array index; it is naively added
   to b, so should be in bytes. */
static gmx_inline gmx_simd_bool_t gmx_load_interaction_mask_pb(long a, const real *b)
{
#ifdef NDEBUG
    return vec_ld(a, (real *) b);
#else
    return vec_lda(a, (real *) b);
#endif
}

#endif /* USE_FUNCTIONS_FOR_QPX */

#endif /* _nbnxn_kernel_simd_utils_ibm_qpx_h_ */
