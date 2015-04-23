/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
#ifndef _nbnxn_kernel_simd_utils_x86_128d_h_
#define _nbnxn_kernel_simd_utils_x86_128d_h_

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/* This files contains all functions/macros for the SIMD kernels
 * which have explicit dependencies on the j-cluster size and/or SIMD-width.
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

#define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

typedef SimdInt32 gmx_exclfilter;
/* This is set to a constant for now, since the code does not adapt automatically just
 * because we set the SIMD widths to other values.
 */
static const int filter_stride = 2;

/* Transpose 2 double precision registers */
static gmx_inline void
gmx_mm_transpose2_op_pd(SimdDouble in0, SimdDouble in1,
                        SimdDouble *out0, SimdDouble *out1)
{
    out0->r = _mm_unpacklo_pd(in0.r, in1.r);
    out1->r = _mm_unpackhi_pd(in0.r, in1.r);
}

/* Sum the elements within each input register and store the sums in out */
static gmx_inline SimdDouble
gmx_mm_transpose_sum2_pr(SimdDouble in0, SimdDouble in1)
{
    SimdDouble tr0, tr1;

    gmx_mm_transpose2_op_pd(in0, in1, &tr0, &tr1);

    return simdAdd(tr0, tr1);
}

static gmx_inline void
load_lj_pair_params(const real *nbfp, const int *type, int aj,
                    SimdDouble *c6_S, SimdDouble *c12_S)
{
    SimdDouble clj_S[UNROLLJ];
    int        p;

    for (p = 0; p < UNROLLJ; p++)
    {
        clj_S[p] = simdLoad(nbfp+type[aj+p]*nbfp_stride);
    }
    gmx_mm_transpose2_op_pd(clj_S[0], clj_S[1], c6_S, c12_S);
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
 * Because of this, we always align the table buffer and provide it in the ti
 * parameter here, even though it is only used with full-width AVX_256. */

static gmx_inline void
load_table_f(const real *tab_coul_F, SimdInt32 ti_S, int gmx_unused *ti,
             SimdDouble *ctab0_S, SimdDouble *ctab1_S)
{
    int        idx[2];
    SimdDouble ctab_S[2];

    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = simdExtractI<0>(ti_S);
    ctab_S[0] = simdLoadU(tab_coul_F+idx[0]);
    idx[1]    = simdExtractI<1>(ti_S);
    ctab_S[1] = simdLoadU(tab_coul_F+idx[1]);

    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], ctab0_S, ctab1_S);
    /* The second force table entry should contain the difference */
    *ctab1_S = simdSub(*ctab1_S, *ctab0_S);
}

static gmx_inline void
load_table_f_v(const real *tab_coul_F, const real *tab_coul_V,
               SimdInt32 ti_S, int gmx_unused *ti,
               SimdDouble *ctab0_S, SimdDouble *ctab1_S, SimdDouble *ctabv_S)
{
    int        idx[2];
    SimdDouble ctab_S[4];

    /* Without SSE4.1 the extract macro needs an immediate: unroll */
    idx[0]    = simdExtractI<0>(ti_S);
    ctab_S[0] = simdLoadU(tab_coul_F+idx[0]);
    idx[1]    = simdExtractI<1>(ti_S);
    ctab_S[1] = simdLoadU(tab_coul_F+idx[1]);

    /* Shuffle the force table entries to a convenient order */
    gmx_mm_transpose2_op_pd(ctab_S[0], ctab_S[1], ctab0_S, ctab1_S);
    /* The second force table entry should contain the difference */
    *ctab1_S = simdSub(*ctab1_S, *ctab0_S);

    ctab_S[2] = simdLoadU(tab_coul_V+idx[0]);
    ctab_S[3] = simdLoadU(tab_coul_V+idx[1]);

    /* Shuffle the energy table entries to a single register */
    ctabv_S->r = _mm_shuffle_pd(ctab_S[2].r, ctab_S[3].r, _MM_SHUFFLE2(0, 0));
}

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return simdSet1I(e);
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    /* For now this has to be an explicit-float load since we use stride==2 */
    SimdFInt32 f = simdLoadFI(reinterpret_cast<const std::int32_t *>(i));
    return {
               f.i
    };
}

static gmx_inline SimdBool
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return {
               _mm_castsi128_pd(_mm_cmpeq_epi32(_mm_andnot_si128(m0.i, m1.i), _mm_setzero_si128()))
    };
}

#endif /* _nbnxn_kernel_simd_utils_x86_s128d_h_ */
