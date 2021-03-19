/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
/*! \file
 *  \brief Define CUDA implementation of nbnxn_gpu_data_mgmt.h
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */
#include "gmxpre.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

// TODO We would like to move this down, but the way NbnxmGpu
//      is currently declared means this has to be before gpu_types.h
#include "nbnxm_cuda_types.h"

// TODO Remove this comment when the above order issue is resolved
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.cuh"
#include "gromacs/gpu_utils/pmalloc.h"
#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/nbnxm_gpu_data_mgmt.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxm_cuda.h"

namespace Nbnxm
{

/* This is a heuristically determined parameter for the Kepler
 * and Maxwell architectures for the minimum size of ci lists by multiplying
 * this constant with the # of multiprocessors on the current device.
 * Since the maximum number of blocks per multiprocessor is 16, the ideal
 * count for small systems is 32 or 48 blocks per multiprocessor. Because
 * there is a bit of fluctuations in the generated block counts, we use
 * a target of 44 instead of the ideal value of 48.
 */
static unsigned int gpu_min_ci_balanced_factor = 44;

void gpu_init_platform_specific(NbnxmGpu* /* nb */)
{
    /* set the kernel type for the current GPU */
    /* pick L1 cache configuration */
    cuda_set_cacheconfig();
}

void gpu_upload_shiftvec(NbnxmGpu* nb, const nbnxn_atomdata_t* nbatom)
{
    NBAtomData*         adat        = nb->atdat;
    const DeviceStream& localStream = *nb->deviceStreams[InteractionLocality::Local];

    /* only if we have a dynamic box */
    if (nbatom->bDynamicBox || !adat->shiftVecUploaded)
    {
        static_assert(sizeof(adat->shiftVec[0]) == sizeof(nbatom->shift_vec[0]),
                      "Sizes of host- and device-side shift vectors should be the same.");
        copyToDeviceBuffer(&adat->shiftVec,
                           reinterpret_cast<const Float3*>(nbatom->shift_vec.data()),
                           0,
                           SHIFTS,
                           localStream,
                           GpuApiCallBehavior::Async,
                           nullptr);
        adat->shiftVecUploaded = true;
    }
}

void gpu_free(NbnxmGpu* nb)
{
    if (nb == nullptr)
    {
        return;
    }

    delete nb->timers;
    sfree(nb->timings);

    NBAtomData* atdat   = nb->atdat;
    NBParamGpu* nbparam = nb->nbparam;

    if (nbparam->elecType == ElecType::EwaldTab || nbparam->elecType == ElecType::EwaldTabTwin)
    {
        destroyParamLookupTable(&nbparam->coulomb_tab, nbparam->coulomb_tab_texobj);
    }

    if (!useLjCombRule(nb->nbparam->vdwType))
    {
        destroyParamLookupTable(&nbparam->nbfp, nbparam->nbfp_texobj);
    }

    if (nbparam->vdwType == VdwType::EwaldGeom || nbparam->vdwType == VdwType::EwaldLB)
    {
        destroyParamLookupTable(&nbparam->nbfp_comb, nbparam->nbfp_comb_texobj);
    }

    freeDeviceBuffer(&atdat->shiftVec);
    freeDeviceBuffer(&atdat->fShift);

    freeDeviceBuffer(&atdat->eLJ);
    freeDeviceBuffer(&atdat->eElec);

    freeDeviceBuffer(&atdat->f);
    freeDeviceBuffer(&atdat->xq);
    if (useLjCombRule(nb->nbparam->vdwType))
    {
        freeDeviceBuffer(&atdat->ljComb);
    }
    else
    {
        freeDeviceBuffer(&atdat->atomTypes);
    }

    /* Free plist */
    auto* plist = nb->plist[InteractionLocality::Local];
    freeDeviceBuffer(&plist->sci);
    freeDeviceBuffer(&plist->cj4);
    freeDeviceBuffer(&plist->imask);
    freeDeviceBuffer(&plist->excl);
    delete plist;
    if (nb->bUseTwoStreams)
    {
        auto* plist_nl = nb->plist[InteractionLocality::NonLocal];
        freeDeviceBuffer(&plist_nl->sci);
        freeDeviceBuffer(&plist_nl->cj4);
        freeDeviceBuffer(&plist_nl->imask);
        freeDeviceBuffer(&plist_nl->excl);
        delete plist_nl;
    }

    /* Free nbst */
    pfree(nb->nbst.eLJ);
    nb->nbst.eLJ = nullptr;

    pfree(nb->nbst.eElec);
    nb->nbst.eElec = nullptr;

    pfree(nb->nbst.fShift);
    nb->nbst.fShift = nullptr;

    delete atdat;
    delete nbparam;
    delete nb;

    if (debug)
    {
        fprintf(debug, "Cleaned up CUDA data structures.\n");
    }
}

int gpu_min_ci_balanced(NbnxmGpu* nb)
{
    return nb != nullptr ? gpu_min_ci_balanced_factor * nb->deviceContext_->deviceInfo().prop.multiProcessorCount
                         : 0;
}

void* gpu_get_xq(NbnxmGpu* nb)
{
    assert(nb);

    return static_cast<void*>(nb->atdat->xq);
}

DeviceBuffer<gmx::RVec> gpu_get_f(NbnxmGpu* nb)
{
    assert(nb);

    return reinterpret_cast<DeviceBuffer<gmx::RVec>>(nb->atdat->f);
}

DeviceBuffer<gmx::RVec> gpu_get_fshift(NbnxmGpu* nb)
{
    assert(nb);

    return reinterpret_cast<DeviceBuffer<gmx::RVec>>(nb->atdat->fShift);
}

/* Initialization for X buffer operations on GPU. */
/* TODO  Remove explicit pinning from host arrays from here and manage in a more natural way*/
void nbnxn_gpu_init_x_to_nbat_x(const Nbnxm::GridSet& gridSet, NbnxmGpu* gpu_nbv)
{
    const DeviceStream& localStream   = *gpu_nbv->deviceStreams[InteractionLocality::Local];
    bool                bDoTime       = gpu_nbv->bDoTime;
    const int           maxNumColumns = gridSet.numColumnsMax();

    reallocateDeviceBuffer(&gpu_nbv->cxy_na,
                           maxNumColumns * gridSet.grids().size(),
                           &gpu_nbv->ncxy_na,
                           &gpu_nbv->ncxy_na_alloc,
                           *gpu_nbv->deviceContext_);
    reallocateDeviceBuffer(&gpu_nbv->cxy_ind,
                           maxNumColumns * gridSet.grids().size(),
                           &gpu_nbv->ncxy_ind,
                           &gpu_nbv->ncxy_ind_alloc,
                           *gpu_nbv->deviceContext_);

    for (unsigned int g = 0; g < gridSet.grids().size(); g++)
    {

        const Nbnxm::Grid& grid = gridSet.grids()[g];

        const int  numColumns      = grid.numColumns();
        const int* atomIndices     = gridSet.atomIndices().data();
        const int  atomIndicesSize = gridSet.atomIndices().size();
        const int* cxy_na          = grid.cxy_na().data();
        const int* cxy_ind         = grid.cxy_ind().data();

        reallocateDeviceBuffer(&gpu_nbv->atomIndices,
                               atomIndicesSize,
                               &gpu_nbv->atomIndicesSize,
                               &gpu_nbv->atomIndicesSize_alloc,
                               *gpu_nbv->deviceContext_);

        if (atomIndicesSize > 0)
        {

            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.openTimingRegion(localStream);
            }

            copyToDeviceBuffer(&gpu_nbv->atomIndices,
                               atomIndices,
                               0,
                               atomIndicesSize,
                               localStream,
                               GpuApiCallBehavior::Async,
                               nullptr);

            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.closeTimingRegion(localStream);
            }
        }

        if (numColumns > 0)
        {
            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.openTimingRegion(localStream);
            }

            int* destPtr = &gpu_nbv->cxy_na[maxNumColumns * g];
            copyToDeviceBuffer(
                    &destPtr, cxy_na, 0, numColumns, localStream, GpuApiCallBehavior::Async, nullptr);

            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.closeTimingRegion(localStream);
            }

            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.openTimingRegion(localStream);
            }

            destPtr = &gpu_nbv->cxy_ind[maxNumColumns * g];
            copyToDeviceBuffer(
                    &destPtr, cxy_ind, 0, numColumns, localStream, GpuApiCallBehavior::Async, nullptr);

            if (bDoTime)
            {
                gpu_nbv->timers->xf[AtomLocality::Local].nb_h2d.closeTimingRegion(localStream);
            }
        }
    }

    if (gpu_nbv->bUseTwoStreams)
    {
        // The above data is transferred on the local stream but is a
        // dependency of the nonlocal stream (specifically the nonlocal X
        // buf ops kernel).  We therefore set a dependency to ensure
        // that the nonlocal stream waits on the local stream here.
        // This call records an event in the local stream:
        gpu_nbv->misc_ops_and_local_H2D_done.markEvent(
                *gpu_nbv->deviceStreams[Nbnxm::InteractionLocality::Local]);
        // ...and this call instructs the nonlocal stream to wait on that event:
        gpu_nbv->misc_ops_and_local_H2D_done.enqueueWaitEvent(
                *gpu_nbv->deviceStreams[Nbnxm::InteractionLocality::NonLocal]);
    }

    return;
}

} // namespace Nbnxm
