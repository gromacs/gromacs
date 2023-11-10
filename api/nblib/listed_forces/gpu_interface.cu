/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Listed forces GPU implementation
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "nblib/listed_forces/gpu_interface.h"

#include <type_traits>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"

#include "nblib/box.h"
#include "nblib/listed_forces/aggregate_transformations.hpp"
#include "nblib/listed_forces/conversionscommon.h"
#include "nblib/listed_forces/definitions.h"
#include "nblib/listed_forces/gpu_kernel.cuh"
#include "nblib/listed_forces/traits.h"
#include "nblib/pbc.hpp"
#include "nblib/util/array.hpp"
#include "nblib/util/traits.hpp"
#include "nblib/util/util.hpp"

namespace gmx
{

ListedForcesNblibGpuImpl::ListedForcesNblibGpuImpl(const gmx_ffparams_t& ffparams,
                                                   float                 electrostaticsScaleFactor,
                                                   const DeviceContext&  deviceContext,
                                                   const DeviceStream&   deviceStream,
                                                   gmx_wallcycle*        wcycle) :
    deviceContext_(deviceContext), deviceStream_(deviceStream), wcycle_(wcycle)
{
    d_potentials_.resize(nblib::GpuListedEnergySize{}, 0);
    h_potentials_.resize(nblib::GpuListedEnergySize{}, 0);
}

ListedForcesNblibGpuImpl::~ListedForcesNblibGpuImpl() = default;

void ListedForcesNblibGpuImpl::updateInteractionListsAndDeviceBuffers(ArrayRef<const int> nbnxnAtomOrder,
                                                                      const InteractionDefinitions& idef,
                                                                      void*              xqDevice,
                                                                      DeviceBuffer<RVec> forceDevice,
                                                                      DeviceBuffer<RVec> fshiftDevice)
{
    nblib::ListedInteractionData interactionsCpu = nblib::convertToNblibInteractions(idef);
    nblib::createAggregates(interactionsCpu);
    interactions_.update(interactionsCpu);

    d_xyzq_        = static_cast<util::array<real, 4>*>(xqDevice);
    d_forces_      = reinterpret_cast<util::array<real, 3>*>(forceDevice);
    d_shiftForces_ = reinterpret_cast<util::array<real, 3>*>(fshiftDevice);

    int numInteractionsTot = interactions_.numInteractionsTot();

    kernelLaunchConfig_.blockSize[0] = nblib::c_threadsPerBlock;
    kernelLaunchConfig_.blockSize[1] = 1;
    kernelLaunchConfig_.blockSize[2] = 1;
    kernelLaunchConfig_.gridSize[0]  = nblib::iceil(numInteractionsTot, nblib::c_threadsPerBlock);
    kernelLaunchConfig_.gridSize[1]  = 1;
    kernelLaunchConfig_.gridSize[2]  = 1;
}

void ListedForcesNblibGpuImpl::setPbc(PbcType pbcType, const matrix box, bool canMoleculeSpanPbc)
{
    pbc_ = nblib::PbcHolderAiuc(pbcType, nblib::Box(box[0][0], box[1][1], box[2][2]));
}

template<bool calcVir, bool calcEner>
void ListedForcesNblibGpuImpl::launchKernel()
{
    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuBonded);

    constexpr int sharedMemSize          = calcVir ? sizeof(real) * 3 * gmx::c_numShiftVectors : 0;
    kernelLaunchConfig_.sharedMemorySize = sharedMemSize;

    using ShiftForceType = std::conditional_t<calcVir, util::array<real, 3>, std::nullptr_t>;
    auto kernelPtr = nblib::computeListedForces<calcEner, real, ShiftForceType, nblib::PbcHolderAiuc>;

    const int* deviceBlockTypes = interactions_.deviceBlockTypes();
    const int* deviceIndices    = interactions_.deviceIndices();

    const void* const* devPA          = interactions_.deviceParametersA();
    const void* const* devPB          = interactions_.deviceParametersB();
    const int*         devNumI        = interactions_.deviceNumInteractions();
    const int*         devOffset      = interactions_.deviceInteractionOffsets();
    const int*         devIndexOffset = interactions_.deviceIndexOffsets();

    real*      potentials = thrust::raw_pointer_cast(d_potentials_.data());
    const auto kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      kernelLaunchConfig_,
                                                      &devNumI,
                                                      &deviceBlockTypes,
                                                      &deviceIndices,
                                                      &devPA,
                                                      &devPB,
                                                      &devOffset,
                                                      &devIndexOffset,
                                                      &d_xyzq_,
                                                      &d_forces_,
                                                      &d_shiftForces_,
                                                      &potentials,
                                                      &pbc_);

    launchGpuKernel(kernelPtr, kernelLaunchConfig_, deviceStream_, nullptr, "computeListedForces", kernelArgs);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuBonded);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

bool ListedForcesNblibGpuImpl::haveInteractions() const
{
    return true;
}

void ListedForcesNblibGpuImpl::launchEnergyTransfer()
{
    h_potentials_ = d_potentials_;
}

void ListedForcesNblibGpuImpl::waitAccumulateEnergyTerms(gmx_enerdata_t* enerd)
{
    auto transferEnergy = [enerd, this](auto interaction) {
        using InteractionType    = std::decay_t<decltype(interaction)>;
        constexpr int nblibIndex = nblib::FindIndex<InteractionType, nblib::GpuListedTypes>{};
        constexpr int gmxIndex   = nblib::FindIndex<InteractionType, nblib::GmxToNblibMapping>{};
        enerd->term[gmxIndex] += this->h_potentials_[nblibIndex];
    };

    using TypeLoop = nblib::Reduce<std::tuple, nblib::GpuListedTypes>;
    nblib::for_each_tuple(transferEnergy, TypeLoop{});

    enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::LJ14][0] +=
            h_potentials_[nblib::GpuVdwIndex{}];
    enerd->grpp.energyGroupPairTerms[NonBondedEnergyTerms::Coulomb14][0] +=
            h_potentials_[nblib::GpuCoulombIndex{}];
}

void ListedForcesNblibGpuImpl::clearEnergies()
{
    thrust::fill(d_potentials_.begin(), d_potentials_.end(), 0);
}

template void ListedForcesNblibGpuImpl::launchKernel<true, true>();
template void ListedForcesNblibGpuImpl::launchKernel<true, false>();
template void ListedForcesNblibGpuImpl::launchKernel<false, true>();
template void ListedForcesNblibGpuImpl::launchKernel<false, false>();

} // namespace gmx
