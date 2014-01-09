/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
#ifndef _nbnxn_kernel_simd_utils_h_
#define _nbnxn_kernel_simd_utils_h_

#include "gromacs/legacyheaders/types/simple.h"

/*! \brief Provides hardware-specific utility routines for the SIMD kernels.
 *
 * Defines all functions, typedefs, constants and macros that have
 * explicit dependencies on the j-cluster size, precision, or SIMD
 * width. This includes handling diagonal, Newton and topology
 * exclusions.
 *
 * The functionality which depends on the j-cluster size is:
 *   LJ-parameter lookup
 *   force table lookup
 *   energy group pair energy storage
 */

#if !defined GMX_NBNXN_SIMD_2XNN && !defined GMX_NBNXN_SIMD_4XN
#error "Must define an NBNxN kernel flavour before including NBNxN kernel utility functions"
#endif

#ifdef GMX_SIMD_REFERENCE_PLAIN_C

/* Align a stack-based thread-local working array. */
static gmx_inline int *
prepare_table_load_buffer(const int gmx_unused *array)
{
    return NULL;
}

#include "nbnxn_kernel_simd_utils_ref.h"

#else /* GMX_SIMD_REFERENCE_PLAIN_C */

#if defined  GMX_X86_SSE2 && !defined __MIC__
/* Include x86 SSE2 compatible SIMD functions */

/* Set the stride for the lookup of the two LJ parameters from their
 * (padded) array. We use the minimum supported SIMD memory alignment.
 */
#if defined GMX_DOUBLE
static const int nbfp_stride = 2;
#else
static const int nbfp_stride = 4;
#endif

/* Align a stack-based thread-local working array. Table loads on
 * full-width AVX_256 use the array, but other implementations do
 * not. */
static gmx_inline int *
prepare_table_load_buffer(const int gmx_unused *array)
{
#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE
    return gmx_simd_align_int(array);
#else
    return NULL;
#endif
}

#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE

/* With full AVX-256 SIMD, half SIMD-width table loads are optimal */
#if GMX_SIMD_WIDTH_HERE == 8
#define TAB_FDV0
#endif

#ifdef GMX_DOUBLE
#include "nbnxn_kernel_simd_utils_x86_256d.h"
#else  /* GMX_DOUBLE */
#include "nbnxn_kernel_simd_utils_x86_256s.h"
#endif /* GMX_DOUBLE */

#else  /* defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE */

/* We use the FDV0 table layout when we can use aligned table loads */
#if GMX_SIMD_WIDTH_HERE == 4
#define TAB_FDV0
#endif

#ifdef GMX_DOUBLE
#include "nbnxn_kernel_simd_utils_x86_128d.h"
#else  /* GMX_DOUBLE */
#include "nbnxn_kernel_simd_utils_x86_128s.h"
#endif /* GMX_DOUBLE */

#endif /* defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE */

#else  /* GMX_X86_SSE2 */

#if GMX_SIMD_WIDTH_HERE > 4
static const int nbfp_stride = 4;
#else
static const int nbfp_stride = GMX_SIMD_WIDTH_HERE;
#endif

/* We use the FDV0 table layout when we can use aligned table loads */
#if GMX_SIMD_WIDTH_HERE == 4
#define TAB_FDV0
#endif

#ifdef GMX_CPU_ACCELERATION_IBM_QPX
#include "nbnxn_kernel_simd_utils_ibm_qpx.h"
#endif /* GMX_CPU_ACCELERATION_IBM_QPX */

#ifdef __MIC__
#include "nbnxn_kernel_simd_utils_x86_mic.h"
#endif

#endif /* GMX_X86_SSE2 */
#endif /* GMX_SIMD_REFERENCE_PLAIN_C */

#if GMX_SIMD_WIDTH_HERE == 4 && !defined GMX_SIMD_REFERENCE_PLAIN_C
#define gmx_mm_pr4    gmx_mm_pr
#define gmx_load_pr4  gmx_load_pr
#define gmx_store_pr4 gmx_store_pr
#define gmx_add_pr4   gmx_add_pr
#endif

#ifndef HAVE_GMX_SUM_SIMD /* should be defined for arch with hardware reduce */
static gmx_inline real
gmx_sum_simd2(gmx_mm_pr x, real* b)
{
    gmx_store_pr(b, x);
    return b[0]+b[1];
}

#if GMX_SIMD_WIDTH_HERE >= 4
static gmx_inline real
gmx_sum_simd4(gmx_mm_pr4 x, real* b)
{
    gmx_store_pr4(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#endif

#if GMX_SIMD_WIDTH_HERE == 2
static gmx_inline real gmx_sum_simd(gmx_mm_pr x, real* b)
{
    gmx_store_pr(b, x);
    return b[0]+b[1];
}
#elif GMX_SIMD_WIDTH_HERE == 4
static gmx_inline real gmx_sum_simd(gmx_mm_pr x, real* b)
{
    gmx_store_pr(b, x);
    return b[0]+b[1]+b[2]+b[3];
}
#elif GMX_SIMD_WIDTH_HERE == 8
static gmx_inline real gmx_sum_simd(gmx_mm_pr x, real* b)
{
    gmx_store_pr(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7];
}
#elif GMX_SIMD_WIDTH_HERE == 16
/* This is getting ridiculous, SIMD horizontal adds would help,
 * but this is not performance critical (only used to reduce energies)
 */
static gmx_inline real gmx_sum_simd(gmx_mm_pr x, real* b)
{
    gmx_store_pr(b, x);
    return b[0]+b[1]+b[2]+b[3]+b[4]+b[5]+b[6]+b[7]+b[8]+b[9]+b[10]+b[11]+b[12]+b[13]+b[14]+b[15];
}
#else
#error "unsupported kernel configuration"
#endif
#endif //HAVE_GMX_SUM_SIMD

#ifdef UNROLLJ
/* Add energy register to possibly multiple terms in the energy array */
static inline void add_ener_grp(gmx_mm_pr e_S, real *v, const int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_pr v_S;

        v_S = gmx_load_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE);
        gmx_store_pr(v+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE, gmx_add_pr(v_S, e_S));
    }
}
#endif

#if defined GMX_NBNXN_SIMD_2XNN && defined UNROLLJ
/* As add_ener_grp, but for two groups of UNROLLJ/2 stored in
 * a single SIMD register.
 */
static inline void
add_ener_grp_halves(gmx_mm_pr e_S, real *v0, real *v1, const int *offset_jj)
{
    gmx_mm_hpr e_S0, e_S1;
    int        jj;

    gmx_pr_to_2hpr(e_S, &e_S0, &e_S1);

    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(&v_S, v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v0+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_S, e_S0));
    }
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(&v_S, v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2);
        gmx_store_hpr(v1+offset_jj[jj]+jj*GMX_SIMD_WIDTH_HERE/2, gmx_add_hpr(v_S, e_S1));
    }
}
#endif

#endif /* _nbnxn_kernel_simd_utils_h_ */
