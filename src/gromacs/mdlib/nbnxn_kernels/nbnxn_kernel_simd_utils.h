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
#ifndef _nbnxn_kernel_simd_utils_h_
#define _nbnxn_kernel_simd_utils_h_

#include "config.h"

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/simd/simd.h"

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

#ifdef GMX_SIMD_REFERENCE

/* Align a stack-based thread-local working array. */
static gmx_inline int *
prepare_table_load_buffer(const int gmx_unused *array)
{
    return NULL;
}

#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ref.h"

#else /* GMX_SIMD_REFERENCE */

#if defined  GMX_TARGET_X86 && !defined GMX_SIMD_X86_MIC
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
 * 256-bit AVX use the array, but other implementations do not.
 */
static gmx_inline int *
prepare_table_load_buffer(int gmx_unused *array)
{
#if GMX_SIMD_REAL_WIDTH >= 8 || (defined GMX_DOUBLE && GMX_SIMD_REAL_WIDTH >= 4)
    return gmx_simd_align_i(array);
#else
    return NULL;
#endif
}

#ifdef GMX_DOUBLE
#if GMX_SIMD_REAL_WIDTH == 2
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_128d.h"
#else
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_256d.h"
#endif
#else /* GMX_DOUBLE */
/* In single precision aligned FDV0 table loads are optimal */
#define TAB_FDV0
#if GMX_SIMD_REAL_WIDTH == 4
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_128s.h"
#else
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_256s.h"
#endif
#endif /* GMX_DOUBLE */

#else  /* GMX_TARGET_X86 && !GMX_SIMD_X86_MIC */

#if GMX_SIMD_REAL_WIDTH > 4
/* For width>4 we use unaligned loads. And thus we can use the minimal stride */
static const int nbfp_stride = 2;
#else
static const int nbfp_stride = GMX_SIMD_REAL_WIDTH;
#endif

/* We use the FDV0 table layout when we can use aligned table loads */
#if GMX_SIMD_REAL_WIDTH == 4
#define TAB_FDV0
#endif

#ifdef GMX_SIMD_IBM_QPX
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_ibm_qpx.h"
#endif /* GMX_SIMD_IBM_QPX */

#ifdef GMX_SIMD_X86_MIC
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils_x86_mic.h"
#endif

#endif /* GMX_TARGET_X86 && !GMX_SIMD_X86_MIC */

#endif /* GMX_SIMD_REFERENCE */

/* If the simd width is 4, but simd4 instructions are not defined,
 * reuse the simd real type and the four instructions we need.
 */
#if GMX_SIMD_REAL_WIDTH == 4 && \
    !((!defined GMX_DOUBLE && defined GMX_SIMD4_HAVE_FLOAT) || \
    (defined GMX_DOUBLE && defined GMX_SIMD4_HAVE_DOUBLE))
#define gmx_simd4_real_t    gmx_simd_real_t
#define gmx_simd4_load_r    gmx_simd_load_r
#define gmx_simd4_store_r   gmx_simd_store_r
#define gmx_simd4_add_r     gmx_simd_add_r
#define gmx_simd4_reduce_r  gmx_simd_reduce_r
#endif

#ifdef UNROLLJ
/* Add energy register to possibly multiple terms in the energy array */
static gmx_inline void add_ener_grp(gmx_simd_real_t e_S, real *v, const int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_simd_real_t v_S;

        v_S = gmx_simd_load_r(v+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH);
        gmx_simd_store_r(v+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH, gmx_simd_add_r(v_S, e_S));
    }
}
#endif

#if defined GMX_NBNXN_SIMD_2XNN && defined UNROLLJ
/* As add_ener_grp, but for two groups of UNROLLJ/2 stored in
 * a single SIMD register.
 */
static gmx_inline void
add_ener_grp_halves(gmx_simd_real_t e_S, real *v0, real *v1, const int *offset_jj)
{
    gmx_mm_hpr e_S0, e_S1;
    int        jj;

    gmx_pr_to_2hpr(e_S, &e_S0, &e_S1);

    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(&v_S, v0+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH/2);
        gmx_store_hpr(v0+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH/2, gmx_add_hpr(v_S, e_S0));
    }
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        gmx_mm_hpr v_S;

        gmx_load_hpr(&v_S, v1+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH/2);
        gmx_store_hpr(v1+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH/2, gmx_add_hpr(v_S, e_S1));
    }
}
#endif

#endif /* _nbnxn_kernel_simd_utils_h_ */
