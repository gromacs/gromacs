/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "gromacs/utility/basedefinitions.h"
#ifdef CALC_COUL_EWALD
#include "gromacs/math/utilities.h"
#endif

#include "config.h"

#ifndef GMX_SIMD_J_UNROLL_SIZE
#error "Need to define GMX_SIMD_J_UNROLL_SIZE before including the 4xn kernel common header file"
#endif

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    (GMX_SIMD_REAL_WIDTH/GMX_SIMD_J_UNROLL_SIZE)

/* The stride of all the atom data arrays is max(UNROLLI,unrollj) */
#if GMX_SIMD_REAL_WIDTH >= UNROLLI
#define STRIDE     (GMX_SIMD_REAL_WIDTH/GMX_SIMD_J_UNROLL_SIZE)
#else
#define STRIDE     (UNROLLI)
#endif

// TODO: Remove when all kernels are in the gmx namespace
using namespace gmx;

#if !defined GMX_NBNXN_SIMD_2XNN && !defined GMX_NBNXN_SIMD_4XN
#error "Must define an NBNxN kernel flavour before including NBNxN kernel utility functions"
#endif

// We use the FDV0 tables for width==4 (when we can load it in one go), or if we don't have any unaligned loads
#if GMX_SIMD_REAL_WIDTH == 4 || !GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL
#define TAB_FDV0
#endif


#ifdef UNROLLJ
/* Add energy register to possibly multiple terms in the energy array */
static gmx_inline void add_ener_grp(SimdReal e_S, real *v, const int *offset_jj)
{
    int jj;

    /* We need to balance the number of store operations with
     * the rapidly increases number of combinations of energy groups.
     * We add to a temporary buffer for 1 i-group vs 2 j-groups.
     */
    for (jj = 0; jj < (UNROLLJ/2); jj++)
    {
        SimdReal v_S;

        v_S = load<SimdReal>(v+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH);
        store(v+offset_jj[jj]+jj*GMX_SIMD_REAL_WIDTH, v_S + e_S);
    }
}
#endif

#if GMX_SIMD_HAVE_INT32_LOGICAL
typedef SimdInt32    SimdBitMask;
#else
typedef SimdReal     SimdBitMask;
#endif

static gmx_inline void gmx_simdcall
gmx_load_simd_4xn_interactions(int                               excl,
                               SimdBitMask gmx_unused            filter_S0,
                               SimdBitMask gmx_unused            filter_S1,
                               SimdBitMask gmx_unused            filter_S2,
                               SimdBitMask gmx_unused            filter_S3,
                               real gmx_unused                  *simd_interaction_array,
                               SimdBool                         *interact_S0,
                               SimdBool                         *interact_S1,
                               SimdBool                         *interact_S2,
                               SimdBool                         *interact_S3)
{
#if GMX_SIMD_HAVE_INT32_LOGICAL
    /* Load integer interaction mask */
    SimdInt32 mask_pr_S(excl);
    *interact_S0  = cvtIB2B(testBits( mask_pr_S & filter_S0 ));
    *interact_S1  = cvtIB2B(testBits( mask_pr_S & filter_S1 ));
    *interact_S2  = cvtIB2B(testBits( mask_pr_S & filter_S2 ));
    *interact_S3  = cvtIB2B(testBits( mask_pr_S & filter_S3 ));
#elif GMX_SIMD_HAVE_LOGICAL
    union
    {
#if GMX_DOUBLE
        std::int64_t i;
#else
        std::int32_t i;
#endif
        real         r;
    } conv;

    conv.i = excl;
    SimdReal      mask_pr_S(conv.r);

    *interact_S0  = testBits( mask_pr_S & filter_S0 );
    *interact_S1  = testBits( mask_pr_S & filter_S1 );
    *interact_S2  = testBits( mask_pr_S & filter_S2 );
    *interact_S3  = testBits( mask_pr_S & filter_S3 );
#else
    // Neither real or integer bitwise logical operations supported.
    // Load masks from memory instead.
    SimdReal      zero = setZero();
    *interact_S0  = ( zero < load<SimdReal>( simd_interaction_array + GMX_SIMD_REAL_WIDTH*((excl >> (0 * UNROLLJ)) & 0xF) ) );
    *interact_S1  = ( zero < load<SimdReal>( simd_interaction_array + GMX_SIMD_REAL_WIDTH*((excl >> (1 * UNROLLJ)) & 0xF) ) );
    *interact_S2  = ( zero < load<SimdReal>( simd_interaction_array + GMX_SIMD_REAL_WIDTH*((excl >> (2 * UNROLLJ)) & 0xF) ) );
    *interact_S3  = ( zero < load<SimdReal>( simd_interaction_array + GMX_SIMD_REAL_WIDTH*((excl >> (3 * UNROLLJ)) & 0xF) ) );
#endif
}

/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 * The combination rule defines, LJ_COMB_GEOM or LJ_COMB_LB are currently
 * set before calling the kernel function. We might want to move that
 * to inside the n-loop and have a different combination rule for different
 * ci's, as no combination rule gives a 50% performance hit for LJ.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

/* Assumes all LJ parameters are identical */
/* #define FIX_LJ_C */
