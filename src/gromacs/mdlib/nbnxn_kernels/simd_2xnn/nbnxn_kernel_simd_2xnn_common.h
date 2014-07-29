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
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#include "../../nbnxn_consts.h"
#ifdef CALC_COUL_EWALD
#include "gromacs/math/utilities.h"
#endif

#ifndef GMX_SIMD_J_UNROLL_SIZE
#error "Need to define GMX_SIMD_J_UNROLL_SIZE before including the 2xnn kernel common header file"
#endif

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    (GMX_SIMD_REAL_WIDTH/GMX_SIMD_J_UNROLL_SIZE)

/* The stride of all the atom data arrays is equal to half the SIMD width */
#define STRIDE     (GMX_SIMD_REAL_WIDTH/GMX_SIMD_J_UNROLL_SIZE)

#include "../nbnxn_kernel_simd_utils.h"

static gmx_inline void gmx_simdcall
gmx_load_simd_2xnn_interactions(int                  excl,
                                gmx_exclfilter       filter_S0,
                                gmx_exclfilter       filter_S2,
                                gmx_simd_bool_t     *interact_S0,
                                gmx_simd_bool_t     *interact_S2)
{
    /* Load integer interaction mask */
    gmx_exclfilter mask_pr_S = gmx_load1_exclfilter(excl);
    *interact_S0  = gmx_checkbitmask_pb(mask_pr_S, filter_S0);
    *interact_S2  = gmx_checkbitmask_pb(mask_pr_S, filter_S2);
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
