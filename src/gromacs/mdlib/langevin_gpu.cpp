/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 *
 * \brief Implements Langevin (SD) integrator using CUDA
 *
 * This file contains implementation of the Langevin (SD) integrator
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "langevin_gpu.h"

#include <assert.h>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"

#include "langevin_gpu_internal.h"

namespace gmx
{

LangevinGpu::LangevinGpu(const DeviceContext& deviceContext,
                         const DeviceStream&  deviceStream,
                         const int            numTempCouplGroups,
                         const float          delta_t,
                         const float*         ref_t,
                         const float*         tau_t) :
    deviceContext_(deviceContext), deviceStream_(deviceStream), numTempCouplGroups_(numTempCouplGroups)
{
    numAtoms_ = 0;

    changePinningPolicy(&sdConstEm_, PinningPolicy::PinnedIfSupported);
    changePinningPolicy(&sdSigmaV_, PinningPolicy::PinnedIfSupported);
    sdConstEm_.resize(numTempCouplGroups);
    sdSigmaV_.resize(numTempCouplGroups);

    for (int i = 0; i < numTempCouplGroups; i++)
    {
        /* Compared to the CPU version we lose precision here by using float instead of double. */
        if (tau_t[i] > 0)
        {
            sdConstEm_[i] = exp(-delta_t / tau_t[i]);
        }
        else
        {
            /* No friction and noise on this group */
            sdConstEm_[i] = 1;
        }

        real kT = gmx::c_boltz * ref_t[i];
        /* The mass is accounted for later, since this differs per atom */
        sdSigmaV_[i] = sqrt(kT * (1 - sdConstEm_[i] * sdConstEm_[i]));
    }
    reallocateDeviceBuffer(
            &d_sdSigmaV_, numTempCouplGroups, &numSdSigmaV_, &numSdSigmaVAlloc_, deviceContext_);
    reallocateDeviceBuffer(
            &d_sdConstEm_, numTempCouplGroups, &numSdConstEm_, &numSdConstEmAlloc_, deviceContext_);

    copyToDeviceBuffer(
            &d_sdSigmaV_, sdSigmaV_.data(), 0, numTempCouplGroups, deviceStream_, GpuApiCallBehavior::Async, nullptr);
    copyToDeviceBuffer(
            &d_sdConstEm_, sdConstEm_.data(), 0, numTempCouplGroups, deviceStream_, GpuApiCallBehavior::Async, nullptr);

    makeDistributionTable();
}

LangevinGpu::~LangevinGpu()
{
    freeDeviceBuffer(&d_inverseMasses_);
    freeDeviceBuffer(&d_tempCouplGroups_);
    freeDeviceBuffer(&d_sdSigmaV_);
    freeDeviceBuffer(&d_sdConstEm_);
    freeDeviceBuffer(&d_distributionTable_);
}

void LangevinGpu::integrate(DeviceBuffer<Float3>       d_x,
                            DeviceBuffer<Float3>       d_xp,
                            DeviceBuffer<Float3>       d_v,
                            const DeviceBuffer<Float3> d_f,
                            const real                 dt,
                            const int                  seed,
                            const int                  step,
                            const SDUpdate             updateType)
{
    GMX_ASSERT(numAtoms_ > 0, "The number of atoms needs to be >0.");
    GMX_ASSERT(updateType == SDUpdate::ForcesOnly || updateType == SDUpdate::FrictionAndNoiseOnly
                       || updateType == SDUpdate::Combined,
               "Unknown SD integrator update type.");

    launchLangevinKernel(numAtoms_,
                         d_x,
                         d_xp,
                         d_v,
                         d_f,
                         d_inverseMasses_,
                         dt,
                         seed,
                         step,
                         d_tempCouplGroups_,
                         d_sdSigmaV_,
                         d_sdConstEm_,
                         d_distributionTable_,
                         updateType,
                         deviceStream_);
}

void LangevinGpu::set(const int                            numAtoms,
                      const ArrayRef<const real>           inverseMasses,
                      const ArrayRef<const unsigned short> tempCouplGroups)
{
    numAtoms_ = numAtoms;

    reallocateDeviceBuffer(
            &d_inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, deviceContext_);
    reallocateDeviceBuffer(
            &d_tempCouplGroups_, numAtoms_, &numTempCouplGroups_, &numTempCouplGroupsAlloc_, deviceContext_);
    reallocateDeviceBuffer(&d_distributionTable_,
                           1 << sc_normalDistributionTableBits,
                           &sizeOfDistributionTable_,
                           &sizeOfDistributionTableAlloc_,
                           deviceContext_);

    copyToDeviceBuffer(
            &d_inverseMasses_, inverseMasses.data(), 0, numAtoms_, deviceStream_, GpuApiCallBehavior::Sync, nullptr);

    if (!tempCouplGroups.empty())
    {
        copyToDeviceBuffer(&d_tempCouplGroups_,
                           tempCouplGroups.data(),
                           0,
                           numAtoms_,
                           deviceStream_,
                           GpuApiCallBehavior::Sync,
                           nullptr);
    }
    else
    {
        changePinningPolicy(&dummyTempCouplGroups_, PinningPolicy::PinnedIfSupported);
        dummyTempCouplGroups_.resize(numAtoms_);
        std::fill(dummyTempCouplGroups_.begin(), dummyTempCouplGroups_.end(), 0);
        copyToDeviceBuffer(&d_tempCouplGroups_,
                           dummyTempCouplGroups_.data(),
                           0,
                           numAtoms_,
                           deviceStream_,
                           GpuApiCallBehavior::Sync,
                           nullptr);
    }

    copyToDeviceBuffer(&d_distributionTable_,
                       distributionTable_.data(),
                       0,
                       1 << sc_normalDistributionTableBits,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
}

void LangevinGpu::makeDistributionTable()
{
    /* Fill the table with the integral of a gaussian distribution, which
     * corresponds to the inverse error function.
     * We avoid integrating a gaussian numerically, since that leads to
     * some loss-of-precision which also accumulates so it is worse for
     * larger indices in the table. */
    constexpr std::size_t tableSize   = 1 << sc_normalDistributionTableBits;
    constexpr std::size_t halfSize    = tableSize / 2;
    constexpr double      invHalfSize = 1.0 / halfSize;

    changePinningPolicy(&distributionTable_, PinningPolicy::PinnedIfSupported);
    distributionTable_.resize(tableSize);
    // Fill in all but the extremal entries of the table
    for (std::size_t i = 0; i < halfSize - 1; i++)
    {
        double r = (i + 0.5) * invHalfSize;
        double x = std::sqrt(2.0) * erfinv(r);

        distributionTable_[halfSize - 1 - i] = -x;
        distributionTable_[halfSize + i]     = x;
    }
    // We want to fill in the extremal table entries with
    // values that make the total variance equal to 1, so
    // measure the variance by summing the squares of the
    // other values of the distribution, starting from the
    // smallest values.
    double sumOfSquares = 0;
    for (std::size_t i = 1; i < halfSize; i++)
    {
        double value = distributionTable_[i];
        sumOfSquares += value * value;
    }
    double missingVariance = 1.0 - 2.0 * sumOfSquares / tableSize;
    GMX_RELEASE_ASSERT(missingVariance > 0,
                       "Incorrect computation of tabulated normal distribution");
    double extremalValue              = std::sqrt(0.5 * missingVariance * tableSize);
    distributionTable_[0]             = -extremalValue;
    distributionTable_[tableSize - 1] = extremalValue;
}

} // namespace gmx
