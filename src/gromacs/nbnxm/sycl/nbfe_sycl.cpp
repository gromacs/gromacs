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
 *  \brief
 *  Launches the SYCL implementation of the nonbonded free-energy kernels
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/utility/template_mp.h"

#include "nbfe_sycl_kernel.h"
#include "nbnxm_sycl_types.h"

namespace gmx
{

// clang-format off
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_8
extern template void launchNbfeKernelHelper<8, false>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeKernelHelper<8, true>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeForeignKernelHelper<8>(NbnxmGpu* nb, InteractionLocality iloc);
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_16
extern template void launchNbfeKernelHelper<16, false>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeKernelHelper<16, true>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeForeignKernelHelper<16>(NbnxmGpu* nb, InteractionLocality iloc);
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_32
extern template void launchNbfeKernelHelper<32, false>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeKernelHelper<32, true>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeForeignKernelHelper<32>(NbnxmGpu* nb, InteractionLocality iloc);
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_64
extern template void launchNbfeKernelHelper<64, false>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeKernelHelper<64, true>(NbnxmGpu* nb, const StepWorkload& stepWork, InteractionLocality iloc);
extern template void launchNbfeForeignKernelHelper<64>(NbnxmGpu* nb, InteractionLocality iloc);
#endif
// clang-format on

template<int subGroupSize>
static void launchNbfeKernels(NbnxmGpu*                      nb,
                              const gmx::SimulationWorkload& simulationWork,
                              const gmx::StepWorkload&       stepWork,
                              const InteractionLocality      iloc)
{
    gmx::dispatchTemplatedFunction(
            [&](auto doCalcEnergies_)
            { launchNbfeKernelHelper<subGroupSize, doCalcEnergies_>(nb, stepWork, iloc); },
            stepWork.computeEnergy);

    if (simulationWork.useGpuForeignNonbondedFE && stepWork.computeDhdl)
    {
        launchNbfeForeignKernelHelper<subGroupSize>(nb, iloc);
    }
}

/*! Launch the nonbonded free energy GPU kernels. */
void gpu_launch_free_energy_kernel(NbnxmGpu*                      nb,
                                   const gmx::SimulationWorkload& simulationWork,
                                   const gmx::StepWorkload&       stepWork,
                                   const InteractionLocality      iloc)
{
    const int subGroupSize = getNbnxmSubGroupSize(nb->deviceContext.deviceInfo(), sc_layoutType);
    switch (subGroupSize)
    {
        // Ensure any changes are in sync with device_management_sycl.cpp and the #if above
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_8
        case 8: launchNbfeKernels<8>(nb, simulationWork, stepWork, iloc); break;
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_16
        case 16: launchNbfeKernels<16>(nb, simulationWork, stepWork, iloc); break;
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_32
        case 32: launchNbfeKernels<32>(nb, simulationWork, stepWork, iloc); break;
#endif
#if SYCL_NBNXM_SUPPORTS_SUBGROUP_SIZE_64
        case 64: launchNbfeKernels<64>(nb, simulationWork, stepWork, iloc); break;
#endif
        default: GMX_RELEASE_ASSERT(false, "Unsupported sub-group size");
    }
}

} // namespace gmx
