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
 *  Data management and kernel launch functions for nbnxm hip.
 *
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/math/functions.h"
#include "gromacs/nbnxm/gpu_common.h"

#include "nbnxm_hip_kernel.h"
#include "nbnxm_hip_kernel_pruneonly.h"
#include "nbnxm_hip_kernel_sci_sort.h"
#include "nbnxm_hip_kernel_sum.h"
#include "nbnxm_hip_types.h"

namespace gmx
{

void gpu_launch_kernel_pruneonly(NbnxmGpu* nb, const InteractionLocality iloc, const int numParts)
{
    std::visit(
            [&](auto&& pairlists)
            {
                auto* plist = pairlists[iloc].get();

                if (plist->haveFreshList)
                {
                    GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

                    /* Set rollingPruningNumParts to signal that it is not set */
                    plist->rollingPruningNumParts = 0;
                }
                else
                {
                    if (plist->rollingPruningNumParts == 0)
                    {
                        plist->rollingPruningNumParts = numParts;
                    }
                    else
                    {
                        GMX_ASSERT(numParts == plist->rollingPruningNumParts,
                                   "It is not allowed to change numParts in between list "
                                   "generation steps");
                    }
                }

                /* Compute the number of list entries to prune in this pass */
                const int numSciInPart = divideRoundUp(plist->numSci, numParts);

                /* Don't launch the kernel if there is no work to do */
                if (numSciInPart <= 0)
                {
                    plist->haveFreshList = false;
                    return;
                }

                launchNbnxmKernelPruneOnly(nb, iloc, &numParts, numSciInPart);
                if (plist->haveFreshList)
                {
                    launchNbnxmKernelSciSort(nb, iloc);
                    plist->haveFreshList                   = false;
                    nb->timers->interaction[iloc].didPrune = true; // Mark that pruning has been done
                }
                else
                {
                    nb->timers->interaction[iloc].didRollingPrune =
                            true; // Mark that rolling pruning has been done
                }

                if (GMX_NATIVE_WINDOWS)
                {
                    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
                    /* Windows: force flushing WDDM queue */
                    hipError_t stat = hipStreamQuery(deviceStream.stream());
                    checkDeviceError(stat,
                                     "hipStreamQuery failed at Windows: force flushing WDDM queue");
                }
            },
            nb->plist);
}

void gpu_launch_kernel(NbnxmGpu* nb, const StepWorkload& stepWork, const InteractionLocality iloc)
{
    const NBParamGpu* nbp = nb->nbparam;
    std::visit(
            [&](auto&& pairlists)
            {
                auto* plist = pairlists[iloc].get();

                if (canSkipNonbondedWork(*nb, iloc))
                {
                    plist->haveFreshList = false;
                    return;
                }

                if (nbp->useDynamicPruning && plist->haveFreshList)
                {
                    // Prunes for rlistOuter and rlistInner, sets plist->haveFreshList=false
                    gpu_launch_kernel_pruneonly(nb, iloc, 1);
                }

                if (plist->numSci == 0)
                {
                    /* Don't launch an empty local kernel */
                    return;
                }


                /* Whether we need to call a combined prune and interaction kernel or just an interaction
                 * kernel. doPrune being true implies we are not using dynamic pruning and are in the first
                 * call to the interaction kernel after a neighbour list step */
                bool doPrune = (plist->haveFreshList && !nb->timers->interaction[iloc].didPrune);

                launchNbnxmKernel(nb, stepWork, iloc, doPrune);
                if (doPrune)
                {
                    launchNbnxmKernelSciSort(nb, iloc);
                }

                if constexpr (sc_useEnergyVirialSeparateDeviceReduction)
                {
                    launchNbnxmKernelSumUp(nb, stepWork, iloc);
                }

                if (GMX_NATIVE_WINDOWS)
                {
                    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
                    /* Windows: force flushing WDDM queue */
                    hipError_t stat = hipStreamQuery(deviceStream.stream());
                    checkDeviceError(stat,
                                     "hipStreamQuery failed at Windows: force flushing WDDM queue");
                }
            },
            nb->plist);
}

} // namespace gmx
