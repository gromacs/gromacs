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
 *  \brief
 *  Data management and kernel launch functions for nbnxm sycl.
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "gromacs/nbnxm/gpu_common.h"
#include "gromacs/utility/exceptions.h"

#include "nbnxm_sycl_kernel.h"
#include "nbnxm_sycl_kernel_pruneonly.h"
#include "nbnxm_sycl_types.h"

namespace Nbnxm
{


void nbnxnInsertNonlocalGpuDependency(NbnxmGpu* nb, const InteractionLocality interactionLocality)
{
    const DeviceStream& deviceStream = *nb->deviceStreams[interactionLocality];
    if (nb->bUseTwoStreams)
    {
        if (interactionLocality == InteractionLocality::Local)
        {
            nb->misc_ops_and_local_H2D_done.markEvent(deviceStream);
        }
        else
        {
            nb->misc_ops_and_local_H2D_done.enqueueWaitEvent(deviceStream);
        }
    }
}


/*! \brief
 * Launch asynchronously the download of nonbonded forces from the GPU
 * (and energies/shift forces if required).
 */
void gpu_launch_cpyback(NbnxmGpu*                nb,
                        struct nbnxn_atomdata_t* nbatom,
                        const gmx::StepWorkload& stepWork,
                        const AtomLocality       atomLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    const InteractionLocality iloc = gpuAtomToInteractionLocality(atomLocality);
    GMX_ASSERT(iloc == InteractionLocality::Local
                       || (iloc == InteractionLocality::NonLocal && nb->bNonLocalStreamDoneMarked == false),
               "Non-local stream is indicating that the copy back event is enqueued at the "
               "beginning of the copy back function.");

    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];
    NBAtomData*         adat         = nb->atdat;

    /* don't launch non-local copy-back if there was no non-local work to do */
    if ((iloc == InteractionLocality::NonLocal) && !haveGpuShortRangeWork(*nb, iloc))
    {
        nb->bNonLocalStreamDoneMarked = false;
        return;
    }

    /* local/nonlocal offset and length used for xq and f */
    auto atomsRange = getGpuAtomRange(adat, atomLocality);

    // With DD the local D2H transfer can only start after the non-local kernel has finished.
    if (iloc == InteractionLocality::Local && nb->bNonLocalStreamDoneMarked)
    {
        nb->nonlocal_done.waitForEvent();
        nb->bNonLocalStreamDoneMarked = false;
    }

    /* DtoH f
     * Skip if buffer ops / reduction is offloaded to the GPU.
     */
    if (!stepWork.useGpuFBufferOps)
    {
        GMX_ASSERT(adat->f.elementSize() == sizeof(Float3),
                   "The size of the force buffer element should be equal to the size of float3.");
        copyFromDeviceBuffer(reinterpret_cast<Float3*>(nbatom->out[0].f.data()) + atomsRange.begin(),
                             &adat->f,
                             atomsRange.begin(),
                             atomsRange.size(),
                             deviceStream,
                             GpuApiCallBehavior::Async,
                             nullptr);
    }

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we want the non-local data
       back first. */
    if (iloc == InteractionLocality::NonLocal)
    {
        nb->nonlocal_done.markEvent(deviceStream);
        nb->bNonLocalStreamDoneMarked = true;
    }

    /* only transfer energies in the local stream */
    if (iloc == InteractionLocality::Local)
    {
        /* DtoH fshift when virial is needed */
        if (stepWork.computeVirial)
        {
            GMX_ASSERT(sizeof(*nb->nbst.fShift) == adat->fShift.elementSize(),
                       "Sizes of host- and device-side shift vector elements should be the same.");
            copyFromDeviceBuffer(
                    nb->nbst.fShift, &adat->fShift, 0, SHIFTS, deviceStream, GpuApiCallBehavior::Async, nullptr);
        }

        /* DtoH energies */
        if (stepWork.computeEnergy)
        {
            GMX_ASSERT(sizeof(*nb->nbst.eLJ) == sizeof(float),
                       "Sizes of host- and device-side LJ energy terms should be the same.");
            copyFromDeviceBuffer(
                    nb->nbst.eLJ, &adat->eLJ, 0, 1, deviceStream, GpuApiCallBehavior::Async, nullptr);
            GMX_ASSERT(sizeof(*nb->nbst.eElec) == sizeof(float),
                       "Sizes of host- and device-side electrostatic energy terms should be the "
                       "same.");
            copyFromDeviceBuffer(
                    nb->nbst.eElec, &adat->eElec, 0, 1, deviceStream, GpuApiCallBehavior::Async, nullptr);
        }
    }
}

void gpu_launch_kernel_pruneonly(NbnxmGpu* nb, const InteractionLocality iloc, const int numParts)
{
    gpu_plist* plist = nb->plist[iloc];

    if (plist->haveFreshList)
    {
        GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

        /* Set rollingPruningNumParts to signal that it is not set */
        plist->rollingPruningNumParts = 0;
        plist->rollingPruningPart     = 0;
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
                       "It is not allowed to change numParts in between list generation steps");
        }
    }

    /* Use a local variable for part and update in plist, so we can return here
     * without duplicating the part increment code.
     */
    const int part = plist->rollingPruningPart;

    plist->rollingPruningPart++;
    if (plist->rollingPruningPart >= plist->rollingPruningNumParts)
    {
        plist->rollingPruningPart = 0;
    }

    /* Compute the number of list entries to prune in this pass */
    const int numSciInPart = (plist->nsci - part) / numParts;

    /* Don't launch the kernel if there is no work to do */
    if (numSciInPart <= 0)
    {
        plist->haveFreshList = false;
        return;
    }

    launchNbnxmKernelPruneOnly(nb, iloc, numParts, part, numSciInPart);

    if (plist->haveFreshList)
    {
        plist->haveFreshList = false;
        nb->didPrune[iloc]   = true; // Mark that pruning has been done
    }
    else
    {
        nb->didRollingPrune[iloc] = true; // Mark that rolling pruning has been done
    }
}

void gpu_launch_kernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const Nbnxm::InteractionLocality iloc)
{
    const NBParamGpu* nbp   = nb->nbparam;
    gpu_plist*        plist = nb->plist[iloc];

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

    if (plist->nsci == 0)
    {
        /* Don't launch an empty local kernel */
        return;
    }

    launchNbnxmKernel(nb, stepWork, iloc);
}

} // namespace Nbnxm
