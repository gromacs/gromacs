/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#ifndef GMX_MDLIB_NBNXN_UTIL_H
#define GMX_MDLIB_NBNXN_UTIL_H

#include <cmath>

#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Returns the base-2 log of n.
 * Generates a fatal error when n is not an integer power of 2.
 */
static gmx_inline int get_2log(int n)
{
    int log2;

    log2 = 0;
    while ((1 << log2) < n)
    {
        log2++;
    }
    if ((1 << log2) != n)
    {
        gmx_fatal(FARGS, "nbnxn na_c (%d) is not a power of 2", n);
    }

    return log2;
}

/* Returns the nbnxn i-cluster size in atoms for the nbnxn kernel type */
static gmx_inline int nbnxn_kernel_to_cluster_i_size(int nb_kernel_type)
{
    switch (nb_kernel_type)
    {
        case nbnxnk4x4_PlainC:
        case nbnxnk4xN_SIMD_4xN:
        case nbnxnk4xN_SIMD_2xNN:
            return NBNXN_CPU_CLUSTER_I_SIZE;
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            /* The cluster size for super/sub lists is only set here.
             * Any value should work for the pair-search and atomdata code.
             * The kernels, of course, might require a particular value.
             */
            return c_nbnxnGpuClusterSize;
        default:
            gmx_incons("unknown kernel type");
    }

    return 0;
}

/* Returns the nbnxn i-cluster size in atoms for the nbnxn kernel type */
static gmx_inline int nbnxn_kernel_to_cluster_j_size(int nb_kernel_type)
{
    int nbnxn_simd_width = 0;
    int cj_size          = 0;

#if GMX_SIMD
    nbnxn_simd_width = GMX_SIMD_REAL_WIDTH;
#endif

    switch (nb_kernel_type)
    {
        case nbnxnk4x4_PlainC:
            cj_size = NBNXN_CPU_CLUSTER_I_SIZE;
            break;
        case nbnxnk4xN_SIMD_4xN:
            cj_size = nbnxn_simd_width;
            break;
        case nbnxnk4xN_SIMD_2xNN:
            cj_size = nbnxn_simd_width/2;
            break;
        case nbnxnk8x8x8_GPU:
        case nbnxnk8x8x8_PlainC:
            cj_size = nbnxn_kernel_to_cluster_i_size(nb_kernel_type);
            break;
        default:
            gmx_incons("unknown kernel type");
    }

    return cj_size;
}


#ifdef __cplusplus
}
#endif

#endif
