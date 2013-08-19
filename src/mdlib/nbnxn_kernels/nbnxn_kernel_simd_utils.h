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
#ifndef _nbnxn_kernel_simd_utils_h_
#define _nbnxn_kernel_simd_utils_h_

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

/* Set the stride for the lookup of the two LJ parameters from their
 * (padded) array.
 * Note that currently only arrays with stride 2 and 4 are available.
 * Since the reference code does not require alignment, we can always use 2.
 */
static const int nbfp_stride = 2;

/* Align a stack-based thread-local working array. */
static gmx_inline int *
prepare_table_load_buffer(const int *array)
{
    return NULL;
}

#include "nbnxn_kernel_simd_utils_ref.h"

#else /* GMX_SIMD_REFERENCE_PLAIN_C */

#ifdef GMX_X86_SSE2
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
prepare_table_load_buffer(const int *array)
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

/*
Berk, 2xnn.c had the following code, but I think it is safe to remove now, given the code immediately above.

#if defined GMX_X86_AVX_256 && !defined GMX_DOUBLE
/ * AVX-256 single precision 2x(4+4) kernel,
 * we can do half SIMD-width aligned FDV0 table loads.
 * /
#define TAB_FDV0
#endif
*/

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

#endif /* GMX_X86_SSE2 */
#endif /* GMX_SIMD_REFERENCE_PLAIN_C */


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
