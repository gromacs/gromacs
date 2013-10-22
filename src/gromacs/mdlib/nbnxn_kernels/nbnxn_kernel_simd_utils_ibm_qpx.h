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
#ifndef _nbnxn_kernel_simd_utils_ibm_qpx_h_
#define _nbnxn_kernel_simd_utils_ibm_qpx_h_

typedef gmx_mm_pr gmx_exclfilter;
static const int filter_stride = 1;

/* The 4xn kernel operates on 4-wide i-force registers */
typedef gmx_mm_pr gmx_mm_pr4;

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

/* Collect all [0123] elements of the 4 inputs to out[0123], respectively */
static gmx_inline void
gmx_transpose_4_ps(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c, gmx_mm_pr d,
                   gmx_mm_pr *out0, gmx_mm_pr *out1,
                   gmx_mm_pr *out2, gmx_mm_pr *out3)
{
    /* Prepare control vectors for swizzling. In its third input,
       vec_perm accepts indices into the effective 8-wide SIMD vector
       created by concatenating its first two inputs. Those indices
       map data from the input vectors to the output vector.

       vec_gpci() converts an octal literal of the indices into the
       correct form for vec_perm() to use. That form is an octal digit
       in bits 0-2 of the mantissa of each double. */
    gmx_mm_pr p6420 = vec_gpci(06420);
    gmx_mm_pr p7531 = vec_gpci(07531);

    /* Four-way swizzle (i.e. transpose) of vectors a = a0a1a2a3, etc. */
    gmx_mm_pr b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_mm_pr b3b1a3a1 = vec_perm(a, b, p7531);
    gmx_mm_pr d2d0c2c0 = vec_perm(c, d, p6420);
    gmx_mm_pr d3d1c3c1 = vec_perm(c, d, p7531);
    *out0 = vec_perm(d2d0c2c0, b2b0a2a0, p7531);
    *out1 = vec_perm(d3d1c3c1, b3b1a3a1, p7531);
    *out2 = vec_perm(d2d0c2c0, b2b0a2a0, p6420);
    *out3 = vec_perm(d3d1c3c1, b3b1a3a1, p6420);
}

/* Collect element 0 and 1 of the 4 inputs to out0 and out1, respectively */
static gmx_inline void
gmx_shuffle_4_ps_fil01_to_2_ps(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c, gmx_mm_pr d,
                               gmx_mm_pr *out0, gmx_mm_pr *out1)
{
    gmx_mm_pr p6420 = vec_gpci(06420);
    gmx_mm_pr p7531 = vec_gpci(07531);

    /* Partial four-way swizzle of vectors a = a0a1a2a3, etc. */
    gmx_mm_pr b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_mm_pr b3b1a3a1 = vec_perm(a, b, p7531);
    gmx_mm_pr d2d0c2c0 = vec_perm(c, d, p6420);
    gmx_mm_pr d3d1c3c1 = vec_perm(c, d, p7531);
    *out0 = vec_perm(d2d0c2c0, b2b0a2a0, p7531);
    *out1 = vec_perm(d3d1c3c1, b3b1a3a1, p7531);
}

/* Collect element 2 of the 4 inputs to out */
static gmx_inline gmx_mm_pr
gmx_shuffle_4_ps_fil2_to_1_ps(gmx_mm_pr a, gmx_mm_pr b, gmx_mm_pr c, gmx_mm_pr d)
{
    gmx_mm_pr p6420 = vec_gpci(06420);

    /* Partial four-way swizzle of vectors a = a0a1a2a3, etc. */
    gmx_mm_pr b2b0a2a0 = vec_perm(a, b, p6420);
    gmx_mm_pr d2d0c2c0 = vec_perm(c, d, p6420);
    return vec_perm(d2d0c2c0, b2b0a2a0, p6420);
}

#ifdef TAB_FDV0
/* Align a stack-based thread-local working array. Table loads on QPX
 * use the array, but most other implementations do not. */
static gmx_inline int *
prepare_table_load_buffer(const int *array)
{
    return gmx_simd_align_int(array);
}

static gmx_inline void
load_table_f(const real *tab_coul_FDV0, gmx_epi32 ti_S, int *ti,
             gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S)
{
#ifdef NDEBUG
    /* Just like 256-bit AVX, we need to use memory to get indices
       into integer registers efficiently. */
    vec_st(ti_S, 0, ti);
#else
    vec_sta(ti_S, 0, ti);
#endif

    /* Here we load 4 aligned reals, but we need just 2 elements of each */
    gmx_mm_pr a = gmx_load_pr(tab_coul_FDV0 + ti[0] * nbfp_stride);
    gmx_mm_pr b = gmx_load_pr(tab_coul_FDV0 + ti[1] * nbfp_stride);
    gmx_mm_pr c = gmx_load_pr(tab_coul_FDV0 + ti[2] * nbfp_stride);
    gmx_mm_pr d = gmx_load_pr(tab_coul_FDV0 + ti[3] * nbfp_stride);

    gmx_shuffle_4_ps_fil01_to_2_ps(a, b, c, d, ctab0_S, ctab1_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_FDV0,
               gmx_epi32 ti_S, int *ti,
               gmx_mm_pr *ctab0_S, gmx_mm_pr *ctab1_S,
               gmx_mm_pr *ctabv_S)
{
#ifdef NDEBUG
    /* Just like 256-bit AVX, we need to use memory to get indices
       into integer registers efficiently. */
    vec_st(ti_S, 0, ti);
#else
    vec_sta(ti_S, 0, ti);
#endif

    /* Here we load 4 aligned reals, but we need just 3 elements of each. */
    gmx_mm_pr a = gmx_load_pr(tab_coul_FDV0 + ti[0] * nbfp_stride);
    gmx_mm_pr b = gmx_load_pr(tab_coul_FDV0 + ti[1] * nbfp_stride);
    gmx_mm_pr c = gmx_load_pr(tab_coul_FDV0 + ti[2] * nbfp_stride);
    gmx_mm_pr d = gmx_load_pr(tab_coul_FDV0 + ti[3] * nbfp_stride);

    gmx_shuffle_4_ps_fil01_to_2_ps(a, b, c, d, ctab0_S, ctab1_S);
    *ctabv_S = gmx_shuffle_4_ps_fil2_to_1_ps(a, b, c, d);
}
#else

/* Not required for BlueGene/Q */

#endif

/* Sum the elements within each input register and store the sums in out.
 */
static gmx_inline gmx_mm_pr
gmx_mm_transpose_sum4_pr(gmx_mm_pr a, gmx_mm_pr b,
                         gmx_mm_pr c, gmx_mm_pr d)
{
    gmx_mm_pr a0b0c0d0, a1b1c1d1, a2b2c2d2, a3b3c3d3;
    gmx_transpose_4_ps(a, b, c, d,
                       &a0b0c0d0,
                       &a1b1c1d1,
                       &a2b2c2d2,
                       &a3b3c3d3);
    /* Now reduce the transposed vectors */
    gmx_mm_pr sum01 = gmx_add_pr(a0b0c0d0, a1b1c1d1);
    gmx_mm_pr sim23 = gmx_add_pr(a2b2c2d2, a3b3c3d3);
    return gmx_add_pr(sum01, sim23);
}

#ifdef GMX_DOUBLE
/* In double precision on x86 it can be faster to first calculate
 * single precision square roots for two double precision registers at
 * once and then use double precision Newton-Raphson iteration to
 * reach full double precision. For QPX, we just wrap the usual
 * reciprocal square roots.
 */
static gmx_inline void
gmx_mm_invsqrt2_pd(gmx_mm_pr in0, gmx_mm_pr in1,
                   gmx_mm_pr *out0, gmx_mm_pr *out1)
{
    *out0 = gmx_invsqrt_pr(in0);
    *out1 = gmx_invsqrt_pr(in1);
}
#endif

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    gmx_mm_pr *c6_S, gmx_mm_pr *c12_S)
{
    /* Here we load 4 aligned reals, but we need just 2 elemnts of each. */
    gmx_mm_pr a = gmx_load_pr(nbfp + type[aj+0] * nbfp_stride);
    gmx_mm_pr b = gmx_load_pr(nbfp + type[aj+1] * nbfp_stride);
    gmx_mm_pr c = gmx_load_pr(nbfp + type[aj+2] * nbfp_stride);
    gmx_mm_pr d = gmx_load_pr(nbfp + type[aj+3] * nbfp_stride);

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
static gmx_inline gmx_mm_pb gmx_load_interaction_mask_pb(long a, const real *b)
{
#ifdef NDEBUG
    return vec_ld(a, (real *) b);
#else
    return vec_lda(a, (real *) b);
#endif
}

#endif /* USE_FUNCTIONS_FOR_QPX */

#endif /* _nbnxn_kernel_simd_utils_ibm_qpx_h_ */
