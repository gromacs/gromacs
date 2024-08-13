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

#include <filesystem>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "nbnxm_enums.h"

namespace gmx
{

/*! \brief Returns the base-2 log of n.
 * *
 * Generates a fatal error when n is not an integer power of 2.
 */
static inline int get_2log(int n)
{
    if (!isPowerOfTwo(n))
    {
        gmx_fatal(FARGS, "nbnxn na_c (%d) is not a power of 2", n);
    }

    return log2I(n);
}

//! The nbnxn i-cluster size in atoms for the given NBNxM kernel type
static constexpr int sc_iClusterSize(const NbnxmKernelType kernelType)
{
    switch (kernelType)
    {
        case NbnxmKernelType::Cpu4x4_PlainC:
        case NbnxmKernelType::Cpu4xN_Simd_4xN:
        case NbnxmKernelType::Cpu4xN_Simd_2xNN: return 4;
        case NbnxmKernelType::Gpu8x8x8:
        case NbnxmKernelType::Cpu8x8x8_PlainC: return c_nbnxnGpuClusterSize;
        case NbnxmKernelType::NotSet:
        case NbnxmKernelType::Count: return 0;
    }

    return 0;
}

/*! \brief The nbnxn j-cluster size in atoms for the given NBNxM kernel type
 *
 * \note When including this file in files compiled for SYCL devices only,
 *       this function can not be called for SIMD kernel types. This is asserted.
 */
static constexpr int sc_jClusterSize(const NbnxmKernelType kernelType)
{
    switch (kernelType)
    {
        case NbnxmKernelType::Cpu4x4_PlainC: return 4;
#if GMX_SIMD
        case NbnxmKernelType::Cpu4xN_Simd_4xN: return GMX_SIMD_REAL_WIDTH;
        case NbnxmKernelType::Cpu4xN_Simd_2xNN: return GMX_SIMD_REAL_WIDTH / 2;
#else
        case NbnxmKernelType::Cpu4xN_Simd_4xN: return 0;
        case NbnxmKernelType::Cpu4xN_Simd_2xNN: return 0;
#endif
        case NbnxmKernelType::Gpu8x8x8:
        case NbnxmKernelType::Cpu8x8x8_PlainC: return c_nbnxnGpuClusterSize / 2;
        case NbnxmKernelType::NotSet:
        case NbnxmKernelType::Count: return 0;
    }

    return 0;
}

/*! \brief Returns whether the pair-list corresponding to nb_kernel_type is simple */
static constexpr bool kernelTypeUsesSimplePairlist(const NbnxmKernelType kernelType)
{
    return (kernelType == NbnxmKernelType::Cpu4x4_PlainC || kernelType == NbnxmKernelType::Cpu4xN_Simd_4xN
            || kernelType == NbnxmKernelType::Cpu4xN_Simd_2xNN);
}

//! Returns whether a SIMD kernel is in use
static constexpr bool kernelTypeIsSimd(const NbnxmKernelType kernelType)
{
    return (kernelType == NbnxmKernelType::Cpu4xN_Simd_4xN
            || kernelType == NbnxmKernelType::Cpu4xN_Simd_2xNN);
}

/*! \brief Returns the increase in pairlist radius when including volume of pairs beyond rlist
 *
 * Due to the cluster size the total volume of the pairlist is (much) more
 * than 4/3*pi*rlist^3. This function returns the increase in radius
 * required to match the volume of the pairlist including the atoms pairs
 * that are beyond rlist.
 *
 * \note This routine does not know which cluster layout is used and assumes the most common one.
 *       Therefore this should only be used to estimates, not for setting a pair list buffer.
 */
real nbnxmPairlistVolumeRadiusIncrease(bool useGpu, real atomDensity);

/*! \brief Returns the effective list radius of the pair-list
 *
 * Due to the cluster size the effective pair-list is longer than
 * that of a simple atom pair-list. This function gives the extra distance.
 *
 * \note This routine does not know which cluster layout is used and assumes the most common one.
 *       Therefore this should only be used to estimates, not for setting a pair list buffer.
 */
real nbnxn_get_rlist_effective_inc(int clusterSize, const RVec& averageClusterBoundingBox);

} // namespace gmx

#endif
