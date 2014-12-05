/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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

#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_simd_utils.h"

static gmx_inline void gmx_simdcall
gmx_load_simd_4xn_interactions(int                        excl,
                               gmx_exclfilter gmx_unused  filter_S0,
                               gmx_exclfilter gmx_unused  filter_S1,
                               gmx_exclfilter gmx_unused  filter_S2,
                               gmx_exclfilter gmx_unused  filter_S3,
                               real gmx_unused           *simd_interaction_array,
                               gmx_simd_bool_t           *interact_S0,
                               gmx_simd_bool_t           *interact_S1,
                               gmx_simd_bool_t           *interact_S2,
                               gmx_simd_bool_t           *interact_S3)
{
#if defined GMX_SIMD_X86_SSE2_OR_HIGHER || defined GMX_SIMD_REFERENCE
    /* Load integer interaction mask */
    gmx_exclfilter mask_pr_S = gmx_load1_exclfilter(excl);
    *interact_S0  = gmx_checkbitmask_pb(mask_pr_S, filter_S0);
    *interact_S1  = gmx_checkbitmask_pb(mask_pr_S, filter_S1);
    *interact_S2  = gmx_checkbitmask_pb(mask_pr_S, filter_S2);
    *interact_S3  = gmx_checkbitmask_pb(mask_pr_S, filter_S3);
#elif defined GMX_SIMD_IBM_QPX
    const int size = GMX_SIMD_REAL_WIDTH * sizeof(real);
    *interact_S0  = gmx_load_interaction_mask_pb(size*((excl >> (0 * UNROLLJ)) & 0xF), simd_interaction_array);
    *interact_S1  = gmx_load_interaction_mask_pb(size*((excl >> (1 * UNROLLJ)) & 0xF), simd_interaction_array);
    *interact_S2  = gmx_load_interaction_mask_pb(size*((excl >> (2 * UNROLLJ)) & 0xF), simd_interaction_array);
    *interact_S3  = gmx_load_interaction_mask_pb(size*((excl >> (3 * UNROLLJ)) & 0xF), simd_interaction_array);
#else
#error "Need implementation of gmx_load_simd_4xn_interactions"
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
