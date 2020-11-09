/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * Implements a bonded force calculator
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/box.h"
#include "nblib/exception.h"
#include "nblib/pbc.hpp"
#include "nblib/listed_forces/calculator.h"
#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/helpers.hpp"

namespace nblib
{

ListedForceCalculator::~ListedForceCalculator() = default;

ListedForceCalculator::ListedForceCalculator(const ListedInteractionData& interactions,
                                             size_t                       bufferSize,
                                             int                          nthr,
                                             const Box&                   box) :
    numThreads(nthr),
    masterForceBuffer_(bufferSize, Vec3{ 0, 0, 0 }),
    pbcHolder_(std::make_unique<PbcHolder>(box))
{
    // split up the work
    threadedInteractions_ = splitListedWork(interactions, bufferSize, numThreads);

    // set up the reduction buffers
    int threadRange = bufferSize / numThreads;
    for (int i = 0; i < numThreads; ++i)
    {
        int rangeStart = i * threadRange;
        int rangeEnd   = rangeStart + threadRange;
        // the last range is a bit bigger due to integer division
        if (i == numThreads - 1)
        {
            rangeEnd = bufferSize;
        }

        threadedForceBuffers_.push_back(std::make_unique<ForceBuffer<Vec3>>(
                masterForceBuffer_.data(), rangeStart, rangeEnd));
    }
}

void ListedForceCalculator::computeForcesAndEnergies(gmx::ArrayRef<const Vec3> x, bool usePbc)
{
    energyBuffer_.fill(0);
    std::vector<std::array<real, std::tuple_size<ListedInteractionData>::value>> energiesPerThread(numThreads);

#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; ++thread)
    {
        if (usePbc)
        {
            energiesPerThread[thread] = reduceListedForces(
                    threadedInteractions_[thread], x, threadedForceBuffers_[thread].get(), *pbcHolder_);
        }
        else
        {
            energiesPerThread[thread] = reduceListedForces(
                    threadedInteractions_[thread], x, threadedForceBuffers_[thread].get(), NoPbc{});
        }
    }

    // reduce energies
    for (int thread = 0; thread < numThreads; ++thread)
    {
        for (int type = 0; type < int(energyBuffer_.size()); ++type)
        {
            energyBuffer_[type] += energiesPerThread[thread][type];
        }
    }

    // reduce forces
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; ++thread)
    {
        auto& thisBuffer = *threadedForceBuffers_[thread];
        // access outliers from other threads
        for (int otherThread = 0; otherThread < numThreads; ++otherThread)
        {
            auto& otherBuffer = *threadedForceBuffers_[otherThread];
            for (const auto& outlier : otherBuffer)
            {
                int index = outlier.first;
                // check whether this outlier falls within the range of <thread>
                if (thisBuffer.inRange(index))
                {
                    auto force = outlier.second;
                    masterForceBuffer_[index] += force;
                }
            }
        }
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates, gmx::ArrayRef<Vec3> forces, bool usePbc)
{
    if (coordinates.size() != forces.size())
    {
        throw InputException("Coordinates array and force buffer size mismatch");
    }

    // check if the force buffers have the same size
    if (masterForceBuffer_.size() != forces.size())
    {
        throw InputException("Input force buffer size mismatch with listed forces buffer");
    }

    // compute forces and fill in local buffers
    computeForcesAndEnergies(coordinates, usePbc);

    // add forces to output force buffers
    for (int pIndex = 0; pIndex < int(forces.size()); pIndex++)
    {
        forces[pIndex] += masterForceBuffer_[pIndex];
    }
}
void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    EnergyType&               energies,
                                    bool                      usePbc)
{
    compute(coordinates, forces, usePbc);

    energies = energyBuffer_;
}

} // namespace nblib
