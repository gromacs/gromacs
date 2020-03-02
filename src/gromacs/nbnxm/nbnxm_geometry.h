/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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

#ifndef GMX_NBNXM_NBNXM_GEOMETRY_H
#define GMX_NBNXM_NBNXM_GEOMETRY_H

#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"

#include "pairlist.h"


/* Returns the base-2 log of n.
 * Generates a fatal error when n is not an integer power of 2.
 */
static inline int get_2log(int n)
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

namespace Nbnxm
{

/* The nbnxn i-cluster size in atoms for each nbnxn kernel type */
static constexpr gmx::EnumerationArray<KernelType, int> IClusterSizePerKernelType = {
    { 0, c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize,
      c_nbnxnGpuClusterSize, c_nbnxnGpuClusterSize }
};

/* The nbnxn j-cluster size in atoms for each nbnxn kernel type */
static constexpr gmx::EnumerationArray<KernelType, int> JClusterSizePerKernelType = {
    { 0, c_nbnxnCpuIClusterSize,
#if GMX_SIMD
      GMX_SIMD_REAL_WIDTH, GMX_SIMD_REAL_WIDTH / 2,
#else
      0, 0,
#endif
      c_nbnxnGpuClusterSize, c_nbnxnGpuClusterSize / 2 }
};

/* Returns whether the pair-list corresponding to nb_kernel_type is simple */
static inline bool kernelTypeUsesSimplePairlist(const KernelType kernelType)
{
    return (kernelType == KernelType::Cpu4x4_PlainC || kernelType == KernelType::Cpu4xN_Simd_4xN
            || kernelType == KernelType::Cpu4xN_Simd_2xNN);
}

static inline bool kernelTypeIsSimd(const KernelType kernelType)
{
    return (kernelType == KernelType::Cpu4xN_Simd_4xN || kernelType == KernelType::Cpu4xN_Simd_2xNN);
}

} // namespace Nbnxm

/* Returns the effective list radius of the pair-list
 *
 * Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 *
 * NOTE: If the i- and j-cluster sizes are identical and you know
 *       the physical dimensions of the clusters, use the next function
 *       for more accurate results
 */
real nbnxn_get_rlist_effective_inc(int jClusterSize, real atomDensity);

/* Returns the effective list radius of the pair-list
 *
 * Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int clusterSize, const gmx::RVec& averageClusterBoundingBox);

#endif
