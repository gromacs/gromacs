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
 * \brief Backend-agnostic functions of Langevin (SD) integrator for GPUs
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

#include "stochasticdynamics_gpu.h"

#include <assert.h>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/functions.h"

namespace gmx
{

namespace
{

//! Size of the table to use for our distribution.
static constexpr std::size_t tableSize = 1 << sc_normalDistributionTableBits;

//! Fill the distribution table with values for the normal distribution.
HostVector<float> makeDistributionTable(const DeviceContext& deviceContext)
{
    HostVector<float> table(
            tableSize, HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported });
    /* Fill the table with the integral of a gaussian distribution, which
     * corresponds to the inverse error function.
     * We avoid integrating a gaussian numerically, since that leads to
     * some loss-of-precision which also accumulates so it is worse for
     * larger indices in the table. */
    static constexpr std::size_t halfSize    = tableSize / 2;
    static constexpr double      invHalfSize = 1.0 / halfSize;

    // Fill in all but the extremal entries of the table
    for (std::size_t i = 0; i < halfSize - 1; i++)
    {
        double r = (i + 0.5) * invHalfSize;
        double x = std::sqrt(2.0) * erfinv(r);

        table[halfSize - 1 - i] = -x;
        table[halfSize + i]     = x;
    }
    // We want to fill in the extremal table entries with
    // values that make the total variance equal to 1, so
    // measure the variance by summing the squares of the
    // other values of the distribution, starting from the
    // smallest values.
    double sumOfSquares = 0;
    for (std::size_t i = 1; i < halfSize; i++)
    {
        double value = table[i];
        sumOfSquares += value * value;
    }
    double missingVariance = 1.0 - 2.0 * sumOfSquares / tableSize;
    GMX_RELEASE_ASSERT(missingVariance > 0,
                       "Incorrect computation of tabulated normal distribution");
    double extremalValue = std::sqrt(0.5 * missingVariance * tableSize);
    table[0]             = -extremalValue;
    table[tableSize - 1] = extremalValue;
    return table;
}

} // namespace

StochasticDynamicsGpu::StochasticDynamicsGpu(const DeviceContext& deviceContext,
                                             const DeviceStream&  deviceStream,
                                             const int            numTempCouplGroups,
                                             const float          delta_t,
                                             const float*         ref_t,
                                             const float*         tau_t) :
    deviceContext_(deviceContext),
    deviceStream_(deviceStream),
    h_globalAtomIndices_(HostAllocationPolicy{ deviceContext, PinningPolicy::PinnedIfSupported })
{
    numAtoms_ = 0;

    HostVector<float> sdConstEm(numTempCouplGroups,
                                HostAllocationPolicy{ deviceContext_, PinningPolicy::PinnedIfSupported });
    HostVector<float> sdSigmaV(numTempCouplGroups,
                               HostAllocationPolicy{ deviceContext_, PinningPolicy::PinnedIfSupported });
    for (int i = 0; i < numTempCouplGroups; i++)
    {
        /* Compared to the CPU version we lose precision here by using float instead of double. */
        if (tau_t[i] > 0)
        {
            sdConstEm[i] = exp(-delta_t / tau_t[i]);
        }
        else
        {
            /* No friction and noise on this group */
            sdConstEm[i] = 1;
        }

        real kT = gmx::c_boltz * ref_t[i];
        /* The mass is accounted for later, since this differs per atom */
        sdSigmaV[i] = sqrt(kT * (1 - sdConstEm[i] * sdConstEm[i]));
    }
    allocateDeviceBuffer(&d_sdSigmaV_, numTempCouplGroups, deviceContext_);
    allocateDeviceBuffer(&d_sdConstEm_, numTempCouplGroups, deviceContext_);
    allocateDeviceBuffer(&d_distributionTable_, 1 << sc_normalDistributionTableBits, deviceContext_);

    auto distributionTable = makeDistributionTable(deviceContext_);

    copyToDeviceBuffer(
            &d_sdSigmaV_, sdSigmaV.data(), 0, numTempCouplGroups, deviceStream_, GpuApiCallBehavior::Async, nullptr);
    copyToDeviceBuffer(
            &d_sdConstEm_, sdConstEm.data(), 0, numTempCouplGroups, deviceStream_, GpuApiCallBehavior::Async, nullptr);
    copyToDeviceBuffer(&d_distributionTable_,
                       distributionTable.data(),
                       0,
                       1 << sc_normalDistributionTableBits,
                       deviceStream_,
                       GpuApiCallBehavior::Async,
                       nullptr);
    deviceStream_.synchronize();
}

StochasticDynamicsGpu::~StochasticDynamicsGpu()
{
    try
    {
        // Wait for all the tasks to complete before freeing the memory. See #4519.
        deviceStream_.synchronize();
        freeDeviceBuffer(&d_inverseMasses_);
        freeDeviceBuffer(&d_tempCouplGroups_);
        freeDeviceBuffer(&d_globalAtomIndices_);
        freeDeviceBuffer(&d_sdSigmaV_);
        freeDeviceBuffer(&d_sdConstEm_);
        freeDeviceBuffer(&d_distributionTable_);
    }
    catch (gmx::InternalError& e)
    {
        fprintf(stderr, "Internal error in destructor of StochasticDynamicsGpu: %s\n", e.what());
    }
}

void StochasticDynamicsGpu::set(const int                            numAtoms,
                                const ArrayRef<const real>           inverseMasses,
                                const ArrayRef<const unsigned short> tempCouplGroups,
                                const ArrayRef<const int>            globalAtomIndices)
{
    numAtoms_ = numAtoms;

    reallocateDeviceBuffer(
            &d_inverseMasses_, numAtoms_, &numInverseMasses_, &numInverseMassesAlloc_, deviceContext_);
    reallocateDeviceBuffer(
            &d_tempCouplGroups_, numAtoms_, &numTempCouplGroups_, &numTempCouplGroupsAlloc_, deviceContext_);

    copyToDeviceBuffer(
            &d_inverseMasses_, inverseMasses.data(), 0, numAtoms_, deviceStream_, GpuApiCallBehavior::Async, nullptr);

    if (!tempCouplGroups.empty())
    {
        copyToDeviceBuffer(&d_tempCouplGroups_,
                           tempCouplGroups.data(),
                           0,
                           numAtoms_,
                           deviceStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
    }
    else
    {
        clearDeviceBufferAsync(&d_tempCouplGroups_, 0, numAtoms_, deviceStream_);
    }

    haveGlobalAtomIndices_ = !globalAtomIndices.empty();
    if (haveGlobalAtomIndices_)
    {
        GMX_RELEASE_ASSERT(globalAtomIndices.ssize() >= numAtoms_,
                           "Need a global atom index for every local atom");
        h_globalAtomIndices_.assign(globalAtomIndices.begin(), globalAtomIndices.begin() + numAtoms_);
        reallocateDeviceBuffer(
                &d_globalAtomIndices_, numAtoms_, &numGlobalAtomIndices_, &numGlobalAtomIndicesAlloc_, deviceContext_);
        copyToDeviceBuffer(&d_globalAtomIndices_,
                           h_globalAtomIndices_.data(),
                           0,
                           numAtoms_,
                           deviceStream_,
                           GpuApiCallBehavior::Async,
                           nullptr);
    }
}


} // namespace gmx
