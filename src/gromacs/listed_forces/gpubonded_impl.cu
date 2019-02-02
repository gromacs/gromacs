/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Implements GPU bonded lists for CUDA
 *
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include "gpubonded_impl.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_vec.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"

struct t_forcerec;

namespace gmx
{

// ---- GpuBonded::Impl

GpuBonded::Impl::Impl(const gmx_ffparams_t &ffparams,
                      void                 *streamPtr)
{
    stream = *static_cast<CommandStream*>(streamPtr);

    allocateDeviceBuffer(&forceparamsDevice, ffparams.numTypes(), nullptr);
    // This could be an async transfer (if the source is pinned), so
    // long as it uses the same stream as the kernels and we are happy
    // to consume additional pinned pages.
    copyToDeviceBuffer(&forceparamsDevice, ffparams.iparams.data(),
                       0, ffparams.numTypes(),
                       stream, GpuApiCallBehavior::Sync, nullptr);
    vtot.resize(F_NRE);
    allocateDeviceBuffer(&vtotDevice, F_NRE, nullptr);
    clearDeviceBufferAsync(&vtotDevice, 0, F_NRE, stream);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        iListsDevice[ftype].nr     = 0;
        iListsDevice[ftype].iatoms = nullptr;
        iListsDevice[ftype].nalloc = 0;
    }
}

GpuBonded::Impl::~Impl()
{
    for (int ftype : ftypesOnGpu)
    {
        if (iListsDevice[ftype].iatoms)
        {
            freeDeviceBuffer(&iListsDevice[ftype].iatoms);
            iListsDevice[ftype].iatoms = nullptr;
        }
    }

    freeDeviceBuffer(&forceparamsDevice);
    freeDeviceBuffer(&vtotDevice);
}

//! Return whether function type \p ftype in \p idef has perturbed interactions
static bool ftypeHasPerturbedEntries(const t_idef  &idef,
                                     int            ftype)
{
    GMX_ASSERT(idef.ilsort == ilsortNO_FE || idef.ilsort == ilsortFE_SORTED,
               "Perturbed interations should be sorted here");

    const t_ilist &ilist = idef.il[ftype];

    return (idef.ilsort != ilsortNO_FE && ilist.nr_nonperturbed != ilist.nr);
}

//! Converts \p src with atom indices in state order to \p dest in nbnxn order
static void convertIlistToNbnxnOrder(const t_ilist       &src,
                                     HostInteractionList *dest,
                                     int                  numAtomsPerInteraction,
                                     ArrayRef<const int>  nbnxnAtomOrder)
{
    GMX_ASSERT(src.size() == 0 || !nbnxnAtomOrder.empty(), "We need the nbnxn atom order");

    dest->iatoms.resize(src.size());

    // TODO use OpenMP to parallelise this loop
    for (int i = 0; i < src.size(); i += 1 + numAtomsPerInteraction)
    {
        dest->iatoms[i] = src.iatoms[i];
        for (int a = 0; a < numAtomsPerInteraction; a++)
        {
            dest->iatoms[i + 1 + a] = nbnxnAtomOrder[src.iatoms[i + 1 + a]];
        }
    }
}

// TODO Consider whether this function should be a factory method that
// makes an object that is the only one capable of the device
// operations needed for the lifetime of an interaction list. This
// would be harder to misuse than GpuBonded, and exchange the problem
// of naming this method for the problem of what to name the
// BondedDeviceInteractionListHandler type.

//! Divides bonded interactions over threads and GPU
void
GpuBonded::Impl::updateInteractionListsAndDeviceBuffers(ArrayRef<const int>  nbnxnAtomOrder,
                                                        const t_idef        &idef,
                                                        void                *xqDevicePtr,
                                                        void                *forceDevicePtr,
                                                        void                *fshiftDevicePtr)
{
    // TODO wallcycle sub start
    haveInteractions_ = false;

    for (int ftype : ftypesOnGpu)
    {
        auto &iList = iLists[ftype];

        /* Perturbation is not implemented in the GPU bonded kernels.
         * But instead of doing all interactions on the CPU, we can
         * still easily handle the types that have no perturbed
         * interactions on the GPU. */
        if (idef.il[ftype].nr > 0 && !ftypeHasPerturbedEntries(idef, ftype))
        {
            haveInteractions_ = true;

            convertIlistToNbnxnOrder(idef.il[ftype],
                                     &iList,
                                     NRAL(ftype), nbnxnAtomOrder);
        }
        else
        {
            iList.iatoms.clear();
        }

        // Update the device if necessary. This can leave some
        // allocations on the device when the host size decreases to
        // zero, which is OK, since we deallocate everything at the
        // end.
        if (iList.size() > 0)
        {
            t_ilist &iListDevice = iListsDevice[ftype];

            reallocateDeviceBuffer(&iListDevice.iatoms, iList.size(), &iListDevice.nr, &iListDevice.nalloc, nullptr);

            copyToDeviceBuffer(&iListDevice.iatoms, iList.iatoms.data(),
                               0, iList.size(),
                               stream, GpuApiCallBehavior::Async, nullptr);
        }
    }

    xqDevice     = static_cast<float4 *>(xqDevicePtr);
    forceDevice  = static_cast<fvec *>(forceDevicePtr);
    fshiftDevice = static_cast<fvec *>(fshiftDevicePtr);
    // TODO wallcycle sub stop
}

bool
GpuBonded::Impl::haveInteractions() const
{
    return haveInteractions_;
}

void
GpuBonded::Impl::launchEnergyTransfer()
{
    // TODO should wrap with ewcLAUNCH_GPU
    GMX_ASSERT(haveInteractions_, "No GPU bonded interactions, so no energies will be computed, so transfer should not be called");

    float *vtot_h   = vtot.data();
    copyFromDeviceBuffer(vtot_h, &vtotDevice,
                         0, F_NRE,
                         stream, GpuApiCallBehavior::Async, nullptr);
}

void
GpuBonded::Impl::accumulateEnergyTerms(gmx_enerdata_t *enerd)
{
    // TODO should wrap with some kind of wait counter, so not all
    // wait goes in to the "Rest" counter
    GMX_ASSERT(haveInteractions_, "No GPU bonded interactions, so no energies will be computed or transferred, so accumulation should not occur");

    cudaError_t stat = cudaStreamSynchronize(stream);
    CU_RET_ERR(stat, "D2H transfer of bonded energies failed");

    for (int ftype : ftypesOnGpu)
    {
        if (ftype != F_LJ14 && ftype != F_COUL14)
        {
            enerd->term[ftype] += vtot[ftype];
        }
    }

    // Note: We do not support energy groups here
    gmx_grppairener_t *grppener = &enerd->grpp;
    GMX_RELEASE_ASSERT(grppener->nener == 1, "No energy group support for bondeds on the GPU");
    grppener->ener[egLJ14][0]   += vtot[F_LJ14];
    grppener->ener[egCOUL14][0] += vtot[F_COUL14];
}

void
GpuBonded::Impl::clearEnergies()
{
    // TODO should wrap with ewcLAUNCH_GPU
    clearDeviceBufferAsync(&vtotDevice, 0, F_NRE, stream);
}

// ---- GpuBonded

GpuBonded::GpuBonded(const gmx_ffparams_t &ffparams,
                     void                 *streamPtr)
    : impl_(new Impl(ffparams, streamPtr))
{
}

GpuBonded::~GpuBonded() = default;

void
GpuBonded::updateInteractionListsAndDeviceBuffers(ArrayRef<const int>  nbnxnAtomOrder,
                                                  const t_idef        &idef,
                                                  void                *xqDevice,
                                                  void                *forceDevice,
                                                  void                *fshiftDevice)
{
    impl_->updateInteractionListsAndDeviceBuffers
        (nbnxnAtomOrder, idef, xqDevice, forceDevice, fshiftDevice);
}

bool
GpuBonded::haveInteractions() const
{
    return impl_->haveInteractions();
}

void
GpuBonded::launchEnergyTransfer()
{
    impl_->launchEnergyTransfer();
}

void
GpuBonded::accumulateEnergyTerms(gmx_enerdata_t *enerd)
{
    impl_->accumulateEnergyTerms(enerd);
}

void
GpuBonded::clearEnergies()
{
    impl_->clearEnergies();
}

}   // namespace gmx
