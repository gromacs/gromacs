/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Implements a bonded force calculator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/listed_forces/calculator.h"

#include <algorithm>

#include "nblib/box.h"
#include "nblib/exception.h"
#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/dataflowpolarization.hpp"
#include "nblib/listed_forces/dataflowrestraints.hpp"
#include "nblib/listed_forces/helpers.hpp"
#include "nblib/listed_forces/positionrestraints.hpp"
#include "nblib/pbc.hpp"

namespace nblib
{

ListedForceCalculator::~ListedForceCalculator() = default;

ListedForceCalculator::ListedForceCalculator(const ListedInteractionData& interactions,
                                             size_t                       bufferSize,
                                             int                          numThreads,
                                             const Box&                   box,
                                             PbcType                      pbcType,
                                             RefCoordScaling              refcoord_scaling) :
    numThreads(numThreads),
    threadedEnergyBuffers_(numThreads),
    threadedForceBuffers_(numThreads),
    threadedShiftForceBuffers_(numThreads),
    threadedVirialsBuffers_(numThreads),
    pbcType_(pbcType),
    box_(box),
    pbcHolder_(std::make_unique<PbcHolderAiuc>(pbcType_, box_)),
    refcoord_scaling_(refcoord_scaling)
{
    // split up the work
    threadedInteractions_ = splitListedWork(interactions, bufferSize, numThreads);

    // set up the reduction buffers
    int threadRange = (bufferSize + numThreads - 1) / numThreads;
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int i = 0; i < numThreads; ++i)
    {
        int rangeStart = i * threadRange;
        int rangeEnd   = rangeStart + threadRange;
        // the last range is a bit bigger due to integer division
        if (i == numThreads - 1)
        {
            rangeEnd = bufferSize;
        }

        threadedForceBuffers_[i]      = ForceBufferProxy<Vec3>(rangeStart, rangeEnd);
        threadedShiftForceBuffers_[i] = std::vector<Vec3>(gmx::c_numShiftVectors);
        threadedVirialsBuffers_[i]    = std::vector<real>(9);
    }
}

ListedForceCalculator::ListedForceCalculator(const ListedInteractionData& interactions,
                                             size_t                       bufferSize,
                                             int                          numThreads,
                                             const Box&                   box) :
    ListedForceCalculator(interactions, bufferSize, numThreads, box, PbcType::Xyz, RefCoordScaling::Default)
{
}

template<class ShiftForce, class Charges, class Virial, class CenterOfMass>
void ListedForceCalculator::computeForcesAndEnergies(gmx::ArrayRef<const Vec3>               x,
                                                     [[maybe_unused]] gmx::ArrayRef<Charges> charges,
                                                     gmx::ArrayRef<Vec3>                     forces,
                                                     [[maybe_unused]] gmx::ArrayRef<ShiftForce> shiftForces,
                                                     [[maybe_unused]] gmx::ArrayRef<Virial> virials,
                                                     [[maybe_unused]] CenterOfMass          com)
{
    if (x.size() != forces.size())
    {
        throw InputException("Coordinates array and force buffer size mismatch");
    }

    constexpr bool haveShiftForces = !std::is_same_v<ShiftForce, std::nullptr_t>;
    if constexpr (haveShiftForces)
    {
        if (shiftForces.size() != gmx::c_numShiftVectors)
        {
            throw InputException("Shift vectors array size mismatch");
        }
    }

    constexpr bool haveVirials      = !std::is_same_v<Virial, std::nullptr_t>;
    constexpr bool haveCenterOfMass = !std::is_same_v<CenterOfMass, std::nullptr_t>;
    // check if charges are specified in the input, will determine if polarization is computed
    constexpr bool haveCharges = std::is_same_v<Charges, const real>;

#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; ++thread)
    {
        std::conditional_t<haveShiftForces, gmx::ArrayRef<Vec3>, gmx::ArrayRef<std::nullptr_t>> shiftForceBuffer;
        if constexpr (haveShiftForces)
        {
            shiftForceBuffer = gmx::ArrayRef<Vec3>(threadedShiftForceBuffers_[thread]);
            std::fill(shiftForceBuffer.begin(), shiftForceBuffer.end(), Vec3{ 0, 0, 0 });
        }

        std::conditional_t<haveVirials, gmx::ArrayRef<real>, gmx::ArrayRef<std::nullptr_t>> virialsBuffer;
        if constexpr (haveVirials)
        {
            virialsBuffer = threadedVirialsBuffers_[thread];
            std::fill(virialsBuffer.begin(), virialsBuffer.end(), 0);
        }

        ForceBufferProxy<Vec3>* threadBuffer = &threadedForceBuffers_[thread];

        // forces in range of this thread are directly written into the output buffer
        threadBuffer->setMasterBuffer(forces);

        // zero out the outliers in the thread buffer
        threadBuffer->clearOutliers();

        threadedEnergyBuffers_[thread] = reduceListedForces(
                threadedInteractions_[thread], x, threadBuffer, shiftForceBuffer, *pbcHolder_);
        if constexpr (haveVirials && haveCenterOfMass)
        {
            threadedEnergyBuffers_[thread] += reduceRestraints(threadedInteractions_[thread],
                                                               x,
                                                               threadBuffer,
                                                               virialsBuffer,
                                                               *pbcHolder_,
                                                               box_,
                                                               pbcType_,
                                                               refcoord_scaling_,
                                                               com);
        }
        if constexpr (haveCharges)
        {
            threadedEnergyBuffers_[thread] += reducePolarization(
                    threadedInteractions_[thread], x, charges, threadBuffer, shiftForceBuffer, *pbcHolder_);
        }
    }

    // reduce shift forces
    // This is a potential candidate for OMP parallelization, but attention should be paid to the
    // relative costs of thread synchronization overhead vs reduction cost in contexts where the
    // number of threads could be large vs where number of threads could be small
    if constexpr (haveShiftForces)
    {
        for (int i = 0; i < gmx::c_numShiftVectors; ++i)
        {
            Vec3 threadSum{ 0, 0, 0 };
            for (int thread = 0; thread < numThreads; ++thread)
            {
                threadSum += (threadedShiftForceBuffers_[thread])[i];
            }
            shiftForces[i] += threadSum;
        }
    }

    // reduce virial
    if constexpr (haveVirials)
    {
        for (int i = 0; i < int(virials.size()); ++i)
        {
            real threadSum = 0;
            for (int thread = 0; thread < numThreads; ++thread)
            {
                threadSum += (threadedVirialsBuffers_[thread])[i];
            }
            virials[i] += threadSum;
        }
    }

    // reduce energies into buffer for thread 0
    for (int thread = 1; thread < numThreads; ++thread)
    {
        threadedEnergyBuffers_[0] += threadedEnergyBuffers_[thread];
    }

    // reduce forces
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; ++thread)
    {
        auto& thisBuffer = threadedForceBuffers_[thread];
        // access outliers from other threads
        for (int otherThread = 0; otherThread < numThreads; ++otherThread)
        {
            auto& otherBuffer = threadedForceBuffers_[otherThread];
            for (const auto& outlier : otherBuffer)
            {
                int index = outlier.first;
                // check whether this outlier falls within the range of <thread>
                if (thisBuffer.inRange(index))
                {
                    auto force = outlier.second;
                    forces[index] += force;
                }
            }
        }
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<const real> charges,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<Vec3>       shiftForces,
                                    gmx::ArrayRef<real>       energies)
{
    computeForcesAndEnergies(
            coordinates, charges, forces, shiftForces, gmx::ArrayRef<std::nullptr_t>{}, true);
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<Vec3>       shiftForces,
                                    gmx::ArrayRef<real>       energies)
{
    computeForcesAndEnergies(coordinates,
                             gmx::ArrayRef<std::nullptr_t>{},
                             forces,
                             shiftForces,
                             gmx::ArrayRef<std::nullptr_t>{},
                             std::nullptr_t{});
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<real>       virials,
                                    gmx::ArrayRef<real>       energies,
                                    Vec3                      com)
{
    computeForcesAndEnergies(
            coordinates, gmx::ArrayRef<std::nullptr_t>{}, forces, gmx::ArrayRef<std::nullptr_t>{}, virials, com);
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<Vec3>       shiftForces,
                                    gmx::ArrayRef<real>       virials,
                                    gmx::ArrayRef<real>       energies,
                                    Vec3                      com)
{
    computeForcesAndEnergies(
            coordinates, gmx::ArrayRef<std::nullptr_t>{}, forces, shiftForces, virials, com);
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<const real> charges,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<Vec3>       shiftForces,
                                    gmx::ArrayRef<real>       virials,
                                    gmx::ArrayRef<real>       energies,
                                    Vec3                      com)
{
    computeForcesAndEnergies(coordinates, charges, forces, shiftForces, virials, com);
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<real>       energies)
{
    computeForcesAndEnergies(coordinates,
                             gmx::ArrayRef<std::nullptr_t>{},
                             forces,
                             gmx::ArrayRef<std::nullptr_t>{},
                             gmx::ArrayRef<std::nullptr_t>{},
                             std::nullptr_t{});
    if (!energies.empty())
    {
        std::copy(threadedEnergyBuffers_[0].begin(), threadedEnergyBuffers_[0].end(), energies.begin());
    }
}

} // namespace nblib
