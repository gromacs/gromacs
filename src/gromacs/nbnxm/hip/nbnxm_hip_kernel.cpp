/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 *  \brief
 *  NBNXM HIP kernels
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_hip_kernel.h"

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_hip_kernel_utils.h"
#include "nbnxm_hip_types.h"

namespace gmx
{

template<bool hasLargeRegisterPool, bool doPruneNBL, bool doCalcEnergies>
void launchNbnxmKernelHelper(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);

// clang-format off
extern template void launchNbnxmKernelHelper<false, false, false>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<false, false, true>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<false, true, true>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<false, true, false>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);

extern template void launchNbnxmKernelHelper<true, false, false>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<true, false, true>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<true, true, true>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
extern template void launchNbnxmKernelHelper<true, true, false>(NbnxmGpu* nb, const StepWorkload&  stepWork, const InteractionLocality iloc);
// clang-format on

void launchNbnxmKernel(NbnxmGpu* nb, const StepWorkload& stepWork, const InteractionLocality iloc, bool doPrune)
{
    const bool hasLargeRegisterPool = targetHasLargeRegisterPool(nb->deviceContext_->deviceInfo());
    const bool doCalcEnergies       = stepWork.computeEnergy;
    dispatchTemplatedFunction(
            [&](auto hasLargeRegisterPool_, auto doCalcEnergies_, auto doPruneNBL_)
            {
                launchNbnxmKernelHelper<hasLargeRegisterPool_, doPruneNBL_, doCalcEnergies_>(
                        nb, stepWork, iloc);
            },
            hasLargeRegisterPool,
            doCalcEnergies,
            doPrune);
}

} // namespace gmx
