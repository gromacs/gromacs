/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief
 * Helper functions and constants for SYCL NBNXM kernels
 *
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H
#define GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H

#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/nbnxm/pairlistparams.h"

namespace Nbnxm
{

#ifndef GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY
//! \brief Default for the prune kernel's jPacked processing concurrency.
#    define GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY 4
#endif

/*! \brief Prune kernel's jPacked processing concurrency.
 *
 *  The \c GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY macro allows compile-time override.
 */
static constexpr int c_syclPruneKernelJPackedConcurrency = GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY;

/* Convenience constants */
/*! \cond */
// cluster size = number of atoms per cluster.
static constexpr int c_clSize = c_nbnxnGpuClusterSize;
// Square of cluster size.
static constexpr int c_clSizeSq = c_clSize * c_clSize;
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

} // namespace Nbnxm

#endif // GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H
