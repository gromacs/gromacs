/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * \brief
 * Helper functions and constants for SYCL NBNXM kernels
 *
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H
#define GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H

#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/nbnxm/pairlistparams.h"

namespace Nbnxm
{

#ifndef GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY
#    define GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY 4
#endif
/*! \brief Macro defining default for the prune kernel's j4 processing concurrency.
 *
 *  The GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY macro allows compile-time override.
 */
static constexpr int c_syclPruneKernelJ4Concurrency = GMX_NBNXN_PRUNE_KERNEL_J4_CONCURRENCY;

/* Convenience constants */
/*! \cond */
// cluster size = number of atoms per cluster.
static constexpr int c_clSize = c_nbnxnGpuClusterSize;
// j-cluster size after split (4 in the current implementation).
static constexpr int c_splitClSize = c_clSize / c_nbnxnGpuClusterpairSplit;
// i-cluster interaction mask for a super-cluster with all c_nbnxnGpuNumClusterPerSupercluster=8 bits set.
static constexpr unsigned superClInteractionMask = ((1U << c_nbnxnGpuNumClusterPerSupercluster) - 1U);

// 1/sqrt(pi), same value as \c M_FLOAT_1_SQRTPI in other NB kernels.
static constexpr float c_OneOverSqrtPi = 0.564189583547756F;
// 1/6, same value as in other NB kernels.
static constexpr float c_oneSixth = 0.16666667F;
// 1/12, same value as in other NB kernels.
static constexpr float c_oneTwelfth = 0.08333333F;
/*! \endcond */

/* The following functions are necessary because on some versions of Intel OpenCL RT, subgroups
 * do not properly work (segfault or create subgroups of size 1) if used in kernels
 * with non-1-dimensional workgroup. */
//! \brief Convert 3D range to 1D
static inline cl::sycl::range<1> flattenRange(cl::sycl::range<3> range3d)
{
    return cl::sycl::range<1>(range3d.size());
}

//! \brief Convert 3D nd_range to 1D
static inline cl::sycl::nd_range<1> flattenNDRange(cl::sycl::nd_range<3> nd_range3d)
{
    return cl::sycl::nd_range<1>(flattenRange(nd_range3d.get_global_range()),
                                 flattenRange(nd_range3d.get_local_range()));
}

//! \brief Convert flattened 1D index to 3D
template<int rangeX, int rangeY>
static inline cl::sycl::id<3> unflattenId(cl::sycl::id<1> id1d)
{
    constexpr unsigned rangeXY = rangeX * rangeY;
    const unsigned     id      = id1d[0];
    const unsigned     z       = id / rangeXY;
    const unsigned     xy      = id % rangeXY;
    return cl::sycl::id<3>(xy % rangeX, xy / rangeX, z);
}

//! \brief Convenience wrapper to do atomic addition to a global buffer
template<cl::sycl::access::mode Mode, class IndexType>
static inline void atomicFetchAdd(DeviceAccessor<float, Mode> acc, const IndexType idx, const float val)
{
    if (cl::sycl::isnormal(val))
    {
        sycl_2020::atomic_ref<float, sycl_2020::memory_order::relaxed, sycl_2020::memory_scope::device, cl::sycl::access::address_space::global_space>
                fout_atomic(acc[idx]);
        fout_atomic.fetch_add(val);
    }
}

static inline float shuffleDown(float var, unsigned int delta, sycl_2020::sub_group sg)
{
    return sg.shuffle_down(var, delta);
}

static inline float shuffleUp(float var, unsigned int delta, sycl_2020::sub_group sg)
{
    return sg.shuffle_up(var, delta);
}

} // namespace Nbnxm

#endif // GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H
