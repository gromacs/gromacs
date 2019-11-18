/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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

/*! \internal \file
 *
 * \brief
 * Declares the ClusterDistanceKernelType enum
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_CLUSTERDISTANCEKERNELTYPE_H
#define GMX_NBNXM_CLUSTERDISTANCEKERNELTYPE_H

#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxassert.h"

#include "atomdata.h"
#include "pairlistparams.h"

//! The types of kernel for calculating the distance between pairs of atom clusters
enum class ClusterDistanceKernelType : int
{
    CpuPlainC,    //!< Plain-C for CPU list
    CpuSimd_4xM,  //!< SIMD for CPU list for j-cluster size matching the SIMD width
    CpuSimd_2xMM, //!< SIMD for CPU list for j-cluster size matching half the SIMD width
    Gpu           //!< For GPU list, can be either plain-C or SIMD
};

//! Return the cluster distance kernel type given the pairlist type and atomdata
static inline ClusterDistanceKernelType getClusterDistanceKernelType(const PairlistType pairlistType,
                                                                     const nbnxn_atomdata_t& atomdata)
{
    if (pairlistType == PairlistType::HierarchicalNxN)
    {
        return ClusterDistanceKernelType::Gpu;
    }
    else if (atomdata.XFormat == nbatXYZ)
    {
        return ClusterDistanceKernelType::CpuPlainC;
    }
    else if (pairlistType == PairlistType::Simple4x2)
    {
#if GMX_SIMD && GMX_SIMD_REAL_WIDTH == 2
        return ClusterDistanceKernelType::CpuSimd_4xM;
#else
        GMX_RELEASE_ASSERT(false, "Expect 2-wide SIMD with 4x2 list and nbat SIMD layout");
#endif
    }
    else if (pairlistType == PairlistType::Simple4x4)
    {
#if GMX_SIMD && GMX_SIMD_REAL_WIDTH == 4
        return ClusterDistanceKernelType::CpuSimd_4xM;
#elif GMX_SIMD && GMX_SIMD_REAL_WIDTH == 8
        return ClusterDistanceKernelType::CpuSimd_2xMM;
#else
        GMX_RELEASE_ASSERT(false,
                           "Expect 4-wide or 8-wide SIMD with 4x4 list and nbat SIMD layout");
#endif
    }
    else
    {
        GMX_ASSERT(pairlistType == PairlistType::Simple4x8, "Unhandled pairlist type");
#if GMX_SIMD && GMX_SIMD_REAL_WIDTH == 8
        return ClusterDistanceKernelType::CpuSimd_4xM;
#elif GMX_SIMD && GMX_SIMD_REAL_WIDTH == 16
        return ClusterDistanceKernelType::CpuSimd_2xMM;
#else
        GMX_RELEASE_ASSERT(false,
                           "Expect 8-wide or 16-wide SIMD with 4x4 list and nbat SIMD layout");
#endif
    }

    GMX_RELEASE_ASSERT(false, "We should have returned before getting here");
    return ClusterDistanceKernelType::CpuPlainC;
}

#endif
