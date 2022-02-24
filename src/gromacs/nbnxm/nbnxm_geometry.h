/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief
 * Declares the geometry-related functionality
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_NBNXM_GEOMETRY_H
#define GMX_NBNXM_NBNXM_GEOMETRY_H

#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"

#include "pairlist.h"


/*! \brief Returns the base-2 log of n.
 * *
 * Generates a fatal error when n is not an integer power of 2.
 */
static inline int get_2log(int n)
{
    if (!gmx::isPowerOfTwo(n))
    {
        gmx_fatal(FARGS, "nbnxn na_c (%d) is not a power of 2", n);
    }

    return gmx::log2I(n);
}

namespace Nbnxm
{

/*! \brief The nbnxn i-cluster size in atoms for each nbnxn kernel type */
static constexpr gmx::EnumerationArray<KernelType, int> IClusterSizePerKernelType = {
    { 0, c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize, c_nbnxnCpuIClusterSize, c_nbnxnGpuClusterSize, c_nbnxnGpuClusterSize }
};

/*! \brief The nbnxn j-cluster size in atoms for each nbnxn kernel type */
static constexpr gmx::EnumerationArray<KernelType, int> JClusterSizePerKernelType = {
    { 0,
      c_nbnxnCpuIClusterSize,
#if GMX_SIMD
      GMX_SIMD_REAL_WIDTH,
      GMX_SIMD_REAL_WIDTH / 2,
#else
      0,
      0,
#endif
      c_nbnxnGpuClusterSize,
      c_nbnxnGpuClusterSize / 2 }
};

/*! \brief Returns whether the pair-list corresponding to nb_kernel_type is simple */
static inline bool kernelTypeUsesSimplePairlist(const KernelType kernelType)
{
    return (kernelType == KernelType::Cpu4x4_PlainC || kernelType == KernelType::Cpu4xN_Simd_4xN
            || kernelType == KernelType::Cpu4xN_Simd_2xNN);
}

//! Returns whether a SIMD kernel is in use
static inline bool kernelTypeIsSimd(const KernelType kernelType)
{
    return (kernelType == KernelType::Cpu4xN_Simd_4xN || kernelType == KernelType::Cpu4xN_Simd_2xNN);
}

} // namespace Nbnxm

/*! \brief Returns the effective list radius of the pair-list
 *
 * Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 *
 * NOTE: If the i- and j-cluster sizes are identical and you know
 *       the physical dimensions of the clusters, use the next function
 *       for more accurate results
 */
real nbnxn_get_rlist_effective_inc(int jClusterSize, real atomDensity);

/*! \brief Returns the effective list radius of the pair-list
 *
 * Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 */
real nbnxn_get_rlist_effective_inc(int clusterSize, const gmx::RVec& averageClusterBoundingBox);

#endif
