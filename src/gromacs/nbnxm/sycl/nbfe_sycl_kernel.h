/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Declares the SYCL nonbonded free-energy kernel launchers
 *
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_H
#define GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_H

#include "nbnxm_sycl_kernel.h"

namespace gmx
{
enum class InteractionLocality;
class StepWorkload;
struct NbnxmGpu;

/*! \brief Number of sub-groups per work-group in the free-energy kernels.
 *
 * Each sub-group handles one entry of the FEP i-atom list, mirroring the
 * one-warp-per-i-atom mapping of the CUDA implementation.
 */
static constexpr int sc_nbfeSubGroupsPerGroup = 2;

/*! \brief Launch the SYCL nonbonded free-energy kernel for the current lambda.
 *
 * \tparam subGroupSize   Sub-group size to compile the kernel for.
 * \tparam doCalcEnergies Whether energies and dV/dlambda are accumulated.
 * \param nb       Non-bonded parameters and buffers.
 * \param stepWork Workload flags for the current step.
 * \param iloc     Interaction locality.
 */
template<int subGroupSize, bool doCalcEnergies>
void launchNbfeKernelHelper(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);

/*! \brief Launch the SYCL nonbonded free-energy kernel for the foreign lambdas.
 *
 * \tparam subGroupSize Sub-group size to compile the kernel for.
 * \param nb    Non-bonded parameters and buffers.
 * \param iloc  Interaction locality.
 */
template<int subGroupSize>
void launchNbfeForeignKernelHelper(NbnxmGpu* nb, InteractionLocality iloc);

} // namespace gmx

#endif // GMX_NBNXM_SYCL_NBFE_SYCL_KERNEL_H
