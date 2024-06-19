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

#include <cstddef>

#include <algorithm>
#include <type_traits>

#include "listed_forces/dataflow.hpp"
#include "listed_forces/helpers.hpp"

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/box.h"
#include "nblib/exception.h"

#include "pbc.hpp"

namespace nblib
{

ListedForceCalculator::~ListedForceCalculator() = default;

ListedForceCalculator::ListedForceCalculator(const ListedInteractionData& interactions,
                                             size_t                       bufferSize,
                                             int                          nthr,
                                             const Box&                   box) :
    numThreads(nthr),
    threadedForceBuffers_(numThreads),
    threadedShiftForceBuffers_(numThreads),
    pbcHolder_(std::make_unique<PbcHolder>(PbcType::Xyz, box))
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
    }
}

template<class ShiftForce>
void ListedForceCalculator::computeForcesAndEnergies(gmx::ArrayRef<const Vec3> x,
                                                     gmx::ArrayRef<Vec3>       forces,
                                                     [[maybe_unused]] gmx::ArrayRef<ShiftForce> shiftForces,
                                                     bool usePbc)
{
    if (x.size() != forces.size())
    {
        throw InputException("Coordinates array and force buffer size mismatch");
    }

    energyBuffer_.fill(0);
    std::vector<std::array<real, std::tuple_size<ListedInteractionData>::value>> energiesPerThread(numThreads);

    constexpr bool haveShiftForces = !std::is_same_v<ShiftForce, std::nullptr_t>;
    if constexpr (haveShiftForces)
    {
        if (shiftForces.size() != gmx::c_numShiftVectors)
        {
            throw InputException("Shift vectors array size mismatch");
        }
    }

#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int thread = 0; thread < numThreads; ++thread)
    {
        std::conditional_t<haveShiftForces, gmx::ArrayRef<Vec3>, gmx::ArrayRef<std::nullptr_t>> shiftForceBuffer;
        if constexpr (haveShiftForces)
        {
            shiftForceBuffer = gmx::ArrayRef<Vec3>(threadedShiftForceBuffers_[thread]);
            std::fill(shiftForceBuffer.begin(), shiftForceBuffer.end(), Vec3{ 0, 0, 0 });
        }

        ForceBufferProxy<Vec3>* threadBuffer = &threadedForceBuffers_[thread];

        // forces in range of this thread are directly written into the output buffer
        threadBuffer->setMainBuffer(forces);

        // zero out the outliers in the thread buffer
        threadBuffer->clearOutliers();

        if (usePbc)
        {
            energiesPerThread[thread] = reduceListedForces(
                    threadedInteractions_[thread], x, threadBuffer, shiftForceBuffer, *pbcHolder_);
        }
        else
        {
            energiesPerThread[thread] = reduceListedForces(
                    threadedInteractions_[thread], x, threadBuffer, shiftForceBuffer, NoPbc{});
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
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<Vec3>       shiftForces,
                                    gmx::ArrayRef<real>       energies,
                                    bool                      usePbc)
{
    computeForcesAndEnergies(coordinates, forces, shiftForces, usePbc);
    if (!energies.empty())
    {
        std::copy(energyBuffer_.begin(), energyBuffer_.end(), energies.begin());
    }
}

void ListedForceCalculator::compute(gmx::ArrayRef<const Vec3> coordinates,
                                    gmx::ArrayRef<Vec3>       forces,
                                    gmx::ArrayRef<real>       energies,
                                    bool                      usePbc)
{
    computeForcesAndEnergies(coordinates, forces, gmx::ArrayRef<std::nullptr_t>{}, usePbc);
    if (!energies.empty())
    {
        std::copy(energyBuffer_.begin(), energyBuffer_.end(), energies.begin());
    }
}

} // namespace nblib
