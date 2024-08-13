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
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/nbnxm/pairlistparams.h"

namespace gmx
{

/*! \brief Prune kernel's jPacked processing concurrency.
 *
 *  The \c GMX_NBNXN_PRUNE_KERNEL_JPACKED_CONCURRENCY macro allows compile-time override.
 */
static constexpr int c_syclPruneKernelJPackedConcurrency = c_pruneKernelJPackedConcurrency;

/*! \endcond */

/*! \brief Explicit uniform load across the warp
 *
 *  Uses the readfirstlane intrinsic to ensure that uniform loads use
 *  scalar registers and subsequent operations on the results generate
 *  scalar instructions.
 *
 *  Note that some ROCm versions' compilers can figure out that the nbnxm
 *  exclusion indices and imasks are uniform and generate the right instructions,
 *  but others (like 4.5 and 5.0.2) do not.
 */
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__) && GMX_GPU_NB_DISABLE_CLUSTER_PAIR_SPLIT
#    define UNIFORM_LOAD_CLUSTER_PAIR_DATA(x) (__builtin_amdgcn_readfirstlane(x))
#else
#    define UNIFORM_LOAD_CLUSTER_PAIR_DATA(x) (x)
#endif

} // namespace gmx

#endif // GMX_NBNXM_SYCL_NBNXN_SYCL_KERNEL_UTILS_H
