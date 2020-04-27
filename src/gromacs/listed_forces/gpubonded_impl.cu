/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include "gpubonded_impl.h"

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/typecasts.cuh"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/forcefieldparameters.h"

struct t_forcerec;

namespace gmx
{
// Number of CUDA threads in a block
constexpr static int c_threadsPerBlock = 256;

// ---- GpuBonded::Impl

GpuBonded::Impl::Impl(const gmx_ffparams_t& ffparams,
                      const float           electrostaticsScaleFactor,
                      const DeviceContext&  deviceContext,
                      const DeviceStream&   deviceStream,
                      gmx_wallcycle*        wcycle) :
    deviceContext_(deviceContext),
    deviceStream_(deviceStream)
{
    GMX_RELEASE_ASSERT(deviceStream.isValid(),
                       "Can't run GPU version of bonded forces in stream that is not valid.");

    static_assert(c_threadsPerBlock >= SHIFTS,
                  "Threads per block in GPU bonded must be >= SHIFTS for the virial kernel "
                  "(calcVir=true)");

    wcycle_ = wcycle;

    allocateDeviceBuffer(&d_forceParams_, ffparams.numTypes(), deviceContext_);
    // This could be an async transfer (if the source is pinned), so
    // long as it uses the same stream as the kernels and we are happy
    // to consume additional pinned pages.
    copyToDeviceBuffer(&d_forceParams_, ffparams.iparams.data(), 0, ffparams.numTypes(),
                       deviceStream_, GpuApiCallBehavior::Sync, nullptr);
    vTot_.resize(F_NRE);
    allocateDeviceBuffer(&d_vTot_, F_NRE, deviceContext_);
    clearDeviceBufferAsync(&d_vTot_, 0, F_NRE, deviceStream_);

    kernelParams_.electrostaticsScaleFactor = electrostaticsScaleFactor;
    kernelParams_.d_forceParams             = d_forceParams_;
    kernelParams_.d_xq                      = d_xq_;
    kernelParams_.d_f                       = d_f_;
    kernelParams_.d_fShift                  = d_fShift_;
    kernelParams_.d_vTot                    = d_vTot_;
    for (int i = 0; i < numFTypesOnGpu; i++)
    {
        kernelParams_.d_iatoms[i]        = nullptr;
        kernelParams_.fTypeRangeStart[i] = 0;
        kernelParams_.fTypeRangeEnd[i]   = -1;
    }

    int fTypeRangeEnd = kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1];

    kernelLaunchConfig_.blockSize[0] = c_threadsPerBlock;
    kernelLaunchConfig_.blockSize[1] = 1;
    kernelLaunchConfig_.blockSize[2] = 1;
    kernelLaunchConfig_.gridSize[0]  = (fTypeRangeEnd + c_threadsPerBlock) / c_threadsPerBlock;
    kernelLaunchConfig_.gridSize[1]  = 1;
    kernelLaunchConfig_.gridSize[2]  = 1;
}

GpuBonded::Impl::~Impl()
{
    for (int fType : fTypesOnGpu)
    {
        if (d_iLists_[fType].iatoms)
        {
            freeDeviceBuffer(&d_iLists_[fType].iatoms);
            d_iLists_[fType].iatoms = nullptr;
        }
    }

    freeDeviceBuffer(&d_forceParams_);
    freeDeviceBuffer(&d_vTot_);
}

//! Return whether function type \p fType in \p idef has perturbed interactions
static bool fTypeHasPerturbedEntries(const InteractionDefinitions& idef, int fType)
{
    GMX_ASSERT(idef.ilsort == ilsortNO_FE || idef.ilsort == ilsortFE_SORTED,
               "Perturbed interations should be sorted here");

    const InteractionList& ilist = idef.il[fType];

    return (idef.ilsort != ilsortNO_FE && idef.numNonperturbedInteractions[fType] != ilist.size());
}

//! Converts \p src with atom indices in state order to \p dest in nbnxn order
static void convertIlistToNbnxnOrder(const InteractionList& src,
                                     HostInteractionList*   dest,
                                     int                    numAtomsPerInteraction,
                                     ArrayRef<const int>    nbnxnAtomOrder)
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

//! Returns \p input rounded up to the closest multiple of \p factor.
static inline int roundUpToFactor(const int input, const int factor)
{
    GMX_ASSERT(factor > 0, "The factor to round up to must be > 0.");

    int remainder = input % factor;

    if (remainder == 0)
    {
        return (input);
    }
    return (input + (factor - remainder));
}

// TODO Consider whether this function should be a factory method that
// makes an object that is the only one capable of the device
// operations needed for the lifetime of an interaction list. This
// would be harder to misuse than GpuBonded, and exchange the problem
// of naming this method for the problem of what to name the
// BondedDeviceInteractionListHandler type.

/*! Divides bonded interactions over threads and GPU.
 *  The bonded interactions are assigned by interaction type to GPU threads. The intereaction
 *  types are assigned in blocks sized as <warp_size>. The beginning and end (thread index) of each
 *  interaction type are stored in kernelParams_. Pointers to the relevant data structures on the
 *  GPU are also stored in kernelParams_.
 *
 * \todo Use DeviceBuffer for the d_xqPtr.
 */
void GpuBonded::Impl::updateInteractionListsAndDeviceBuffers(ArrayRef<const int> nbnxnAtomOrder,
                                                             const InteractionDefinitions& idef,
                                                             void*                         d_xqPtr,
                                                             DeviceBuffer<RVec>            d_fPtr,
                                                             DeviceBuffer<RVec> d_fShiftPtr)
{
    // TODO wallcycle sub start
    haveInteractions_ = false;
    int fTypesCounter = 0;

    for (int fType : fTypesOnGpu)
    {
        auto& iList = iLists_[fType];

        /* Perturbation is not implemented in the GPU bonded kernels.
         * But instead of doing all interactions on the CPU, we can
         * still easily handle the types that have no perturbed
         * interactions on the GPU. */
        if (!idef.il[fType].empty() && !fTypeHasPerturbedEntries(idef, fType))
        {
            haveInteractions_ = true;

            convertIlistToNbnxnOrder(idef.il[fType], &iList, NRAL(fType), nbnxnAtomOrder);
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
            t_ilist& d_iList = d_iLists_[fType];

            reallocateDeviceBuffer(&d_iList.iatoms, iList.size(), &d_iList.nr, &d_iList.nalloc,
                                   deviceContext_);

            copyToDeviceBuffer(&d_iList.iatoms, iList.iatoms.data(), 0, iList.size(), deviceStream_,
                               GpuApiCallBehavior::Async, nullptr);
        }
        kernelParams_.fTypesOnGpu[fTypesCounter]    = fType;
        kernelParams_.numFTypeIAtoms[fTypesCounter] = iList.size();
        int numBonds = iList.size() / (interaction_function[fType].nratoms + 1);
        kernelParams_.numFTypeBonds[fTypesCounter] = numBonds;
        kernelParams_.d_iatoms[fTypesCounter]      = d_iLists_[fType].iatoms;
        if (fTypesCounter == 0)
        {
            kernelParams_.fTypeRangeStart[fTypesCounter] = 0;
        }
        else
        {
            kernelParams_.fTypeRangeStart[fTypesCounter] =
                    kernelParams_.fTypeRangeEnd[fTypesCounter - 1] + 1;
        }
        kernelParams_.fTypeRangeEnd[fTypesCounter] = kernelParams_.fTypeRangeStart[fTypesCounter]
                                                     + roundUpToFactor(numBonds, warp_size) - 1;

        GMX_ASSERT(numBonds > 0
                           || kernelParams_.fTypeRangeEnd[fTypesCounter]
                                      <= kernelParams_.fTypeRangeStart[fTypesCounter],
                   "Invalid GPU listed forces setup. numBonds must be > 0 if there are threads "
                   "allocated to do work on that interaction function type.");
        GMX_ASSERT(kernelParams_.fTypeRangeStart[fTypesCounter] % warp_size == 0
                           && (kernelParams_.fTypeRangeEnd[fTypesCounter] + 1) % warp_size == 0,
                   "The bonded interactions must be assigned to the GPU in blocks of warp size.");

        fTypesCounter++;
    }

    int fTypeRangeEnd               = kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1];
    kernelLaunchConfig_.gridSize[0] = (fTypeRangeEnd + c_threadsPerBlock) / c_threadsPerBlock;

    d_xq_     = static_cast<float4*>(d_xqPtr);
    d_f_      = asFloat3(d_fPtr);
    d_fShift_ = asFloat3(d_fShiftPtr);

    kernelParams_.d_xq          = d_xq_;
    kernelParams_.d_f           = d_f_;
    kernelParams_.d_fShift      = d_fShift_;
    kernelParams_.d_forceParams = d_forceParams_;
    kernelParams_.d_vTot        = d_vTot_;

    // TODO wallcycle sub stop
}

void GpuBonded::Impl::setPbc(PbcType pbcType, const matrix box, bool canMoleculeSpanPbc)
{
    PbcAiuc pbcAiuc;
    setPbcAiuc(canMoleculeSpanPbc ? numPbcDimensions(pbcType) : 0, box, &pbcAiuc);
    kernelParams_.pbcAiuc = pbcAiuc;
}

bool GpuBonded::Impl::haveInteractions() const
{
    return haveInteractions_;
}

void GpuBonded::Impl::launchEnergyTransfer()
{
    GMX_ASSERT(haveInteractions_,
               "No GPU bonded interactions, so no energies will be computed, so transfer should "
               "not be called");
    wallcycle_sub_start_nocount(wcycle_, ewcsLAUNCH_GPU_BONDED);
    // TODO add conditional on whether there has been any compute (and make sure host buffer doesn't contain garbage)
    float* h_vTot = vTot_.data();
    copyFromDeviceBuffer(h_vTot, &d_vTot_, 0, F_NRE, deviceStream_, GpuApiCallBehavior::Async, nullptr);
    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_GPU_BONDED);
}

void GpuBonded::Impl::waitAccumulateEnergyTerms(gmx_enerdata_t* enerd)
{
    GMX_ASSERT(haveInteractions_,
               "No GPU bonded interactions, so no energies will be computed or transferred, so "
               "accumulation should not occur");

    wallcycle_start(wcycle_, ewcWAIT_GPU_BONDED);
    cudaError_t stat = cudaStreamSynchronize(deviceStream_.stream());
    CU_RET_ERR(stat, "D2H transfer of bonded energies failed");
    wallcycle_stop(wcycle_, ewcWAIT_GPU_BONDED);

    for (int fType : fTypesOnGpu)
    {
        if (fType != F_LJ14 && fType != F_COUL14)
        {
            enerd->term[fType] += vTot_[fType];
        }
    }

    // Note: We do not support energy groups here
    gmx_grppairener_t* grppener = &enerd->grpp;
    GMX_RELEASE_ASSERT(grppener->nener == 1, "No energy group support for bondeds on the GPU");
    grppener->ener[egLJ14][0] += vTot_[F_LJ14];
    grppener->ener[egCOUL14][0] += vTot_[F_COUL14];
}

void GpuBonded::Impl::clearEnergies()
{
    wallcycle_start_nocount(wcycle_, ewcLAUNCH_GPU);
    wallcycle_sub_start_nocount(wcycle_, ewcsLAUNCH_GPU_BONDED);
    clearDeviceBufferAsync(&d_vTot_, 0, F_NRE, deviceStream_);
    wallcycle_sub_stop(wcycle_, ewcsLAUNCH_GPU_BONDED);
    wallcycle_stop(wcycle_, ewcLAUNCH_GPU);
}

// ---- GpuBonded

GpuBonded::GpuBonded(const gmx_ffparams_t& ffparams,
                     const float           electrostaticsScaleFactor,
                     const DeviceContext&  deviceContext,
                     const DeviceStream&   deviceStream,
                     gmx_wallcycle*        wcycle) :
    impl_(new Impl(ffparams, electrostaticsScaleFactor, deviceContext, deviceStream, wcycle))
{
}

GpuBonded::~GpuBonded() = default;

void GpuBonded::updateInteractionListsAndDeviceBuffers(ArrayRef<const int>           nbnxnAtomOrder,
                                                       const InteractionDefinitions& idef,
                                                       void*                         d_xq,
                                                       DeviceBuffer<RVec>            d_f,
                                                       DeviceBuffer<RVec>            d_fShift)
{
    impl_->updateInteractionListsAndDeviceBuffers(nbnxnAtomOrder, idef, d_xq, d_f, d_fShift);
}

void GpuBonded::setPbc(PbcType pbcType, const matrix box, bool canMoleculeSpanPbc)
{
    impl_->setPbc(pbcType, box, canMoleculeSpanPbc);
}

bool GpuBonded::haveInteractions() const
{
    return impl_->haveInteractions();
}

void GpuBonded::setPbcAndlaunchKernel(PbcType                  pbcType,
                                      const matrix             box,
                                      bool                     canMoleculeSpanPbc,
                                      const gmx::StepWorkload& stepWork)
{
    setPbc(pbcType, box, canMoleculeSpanPbc);
    launchKernel(stepWork);
}

void GpuBonded::launchEnergyTransfer()
{
    impl_->launchEnergyTransfer();
}

void GpuBonded::waitAccumulateEnergyTerms(gmx_enerdata_t* enerd)
{
    impl_->waitAccumulateEnergyTerms(enerd);
}

void GpuBonded::clearEnergies()
{
    impl_->clearEnergies();
}

} // namespace gmx
