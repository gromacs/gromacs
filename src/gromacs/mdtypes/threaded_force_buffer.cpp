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
 *
 * \brief This file defines the implementation of ThreadForceBuffer and ThreadedForceBuffer.
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "threaded_force_buffer.h"

#include <cstdio>

#include <algorithm>
#include <array>
#include <string>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief The max thread number is arbitrary, we used a fixed number
 * to avoid memory management.  Using more than 16 threads is probably
 * never useful performance wise. */
static constexpr int s_maxNumThreadsForReduction = 256;

template<typename ForceBufferElementType>
ThreadForceBuffer<ForceBufferElementType>::ThreadForceBuffer(const int  threadIndex,
                                                             const bool useEnergyTerms,
                                                             const int  numEnergyGroups) :
    threadIndex_(threadIndex), shiftForces_(c_numShiftVectors), groupPairEnergies_(numEnergyGroups)
{
    if (useEnergyTerms)
    {
        energyTerms_.resize(F_NRE);
    }
}

template<typename ForceBufferElementType>
void ThreadForceBuffer<ForceBufferElementType>::clearForcesAndEnergies()
{
    constexpr int c_numComponents = sizeof(ForceBufferElementType) / sizeof(real);

    for (int atomIndex : usedBlockIndices_)
    {
        const int bufferIndexBegin = atomIndex * s_reductionBlockSize * c_numComponents;
        const int bufferIndexEnd   = bufferIndexBegin + s_reductionBlockSize * c_numComponents;
        std::fill(forceBuffer_.begin() + bufferIndexBegin, forceBuffer_.begin() + bufferIndexEnd, 0.0_real);
    }

    const RVec zeroVec = { 0.0_real, 0.0_real, 0.0_real };
    std::fill(shiftForces_.begin(), shiftForces_.end(), zeroVec);

    std::fill(energyTerms_.begin(), energyTerms_.end(), 0.0_real);

    for (int i = 0; i < static_cast<int>(NonBondedEnergyTerms::Count); i++)
    {
        for (int j = 0; j < groupPairEnergies_.nener; j++)
        {
            groupPairEnergies_.energyGroupPairTerms[i][j] = 0.0_real;
        }
    }
    for (auto i : keysOf(dvdl_))
    {
        dvdl_[i] = 0;
    }
}

template<typename ForceBufferElementType>
void ThreadForceBuffer<ForceBufferElementType>::resizeBufferAndClearMask(const int numAtoms)
{
    numAtoms_ = numAtoms;

    const int numBlocks = (numAtoms + s_reductionBlockSize - 1) >> s_numReductionBlockBits;

    reductionMask_.resize(numBlocks);

    constexpr size_t c_numComponentsInElement = sizeof(ForceBufferElementType) / sizeof(real);
    int              newNumElements           = numBlocks * s_reductionBlockSize;
    if (c_numComponentsInElement != 4 && newNumElements == numAtoms)
    {
        // Pad with one element to allow 4-wide SIMD loads and stores.
        // Note that actually only one real is needed, but we need a whole element for the ArrayRef.
        newNumElements += 1;
    }
    forceBuffer_.resize(newNumElements * c_numComponentsInElement);

    for (gmx_bitmask_t& mask : reductionMask_)
    {
        bitmask_clear(&mask);
    }
}

template<typename ForceBufferElementType>
void ThreadForceBuffer<ForceBufferElementType>::processMask()
{
    // Now we are done setting the masks, generate the new list of used blocks
    usedBlockIndices_.clear();
    for (int b = 0; b < gmx::ssize(reductionMask_); b++)
    {
        if (bitmask_is_set(reductionMask_[b], threadIndex_))
        {
            usedBlockIndices_.push_back(b);
        }
    }
}

namespace
{

//! \brief Reduce thread-local force buffers into \p force (does not reduce shift forces)
template<typename ForceBufferElementType>
void reduceThreadForceBuffers(ArrayRef<gmx::RVec> force,
                              ArrayRef<std::unique_ptr<ThreadForceBuffer<ForceBufferElementType>>> threadForceBuffers,
                              ArrayRef<const gmx_bitmask_t> masks,
                              ArrayRef<const int>           usedBlockIndices)
{
    const int numBuffers = threadForceBuffers.size();
    GMX_ASSERT(numBuffers <= s_maxNumThreadsForReduction,
               "There is a limit on the number of buffers we can use for reduction");

    const int numAtoms = threadForceBuffers[0]->size();

    rvec* gmx_restrict f = as_rvec_array(force.data());

    /* This reduction can run on any number of threads, independently of the number of buffers.
     * But if the number of threads matches the number of buffers, which it currently does,
     * the uniform distribution of the touched blocks over nthreads will
     * match the distribution of bonded over threads well in most cases,
     * which means that threads mostly reduce their own data which increases
     * the number of cache hits.
     * Additionally, we should always use the same number of threads in parallel
     * regions in OpenMP, otherwise the performance will degrade significantly.
     */
    const int gmx_unused numThreadsForReduction = threadForceBuffers.size();
// nvc++ 24.1+ version has bug due to which it generates incorrect OMP code for this region
// so disable this until nvc++ gets fixed.
#if !defined(__NVCOMPILER)
#    pragma omp parallel for num_threads(numThreadsForReduction) schedule(static)
#endif
    for (int b = 0; b < usedBlockIndices.ssize(); b++)
    {
        try
        {
            // Reduce the buffers that contribute to this block
            ForceBufferElementType* fp[s_maxNumThreadsForReduction];

            const int blockIndex = usedBlockIndices[b];

            // Make a list of threads that have this block index set in the mask
            int numContributingBuffers = 0;
            for (int ft = 0; ft < numBuffers; ft++)
            {
                if (bitmask_is_set(masks[blockIndex], ft))
                {
                    fp[numContributingBuffers++] =
                            threadForceBuffers[ft]->forceBufferWithPadding().paddedArrayRef().data();
                }
            }
            if (numContributingBuffers > 0)
            {
                // Reduce the selected buffers
                int a0 = blockIndex * ThreadForceBuffer<ForceBufferElementType>::s_reductionBlockSize;
                int a1 = (blockIndex + 1) * ThreadForceBuffer<ForceBufferElementType>::s_reductionBlockSize;
                // Note: It would be nice if we could pad f to avoid this min()
                a1 = std::min(a1, numAtoms);
                if (numContributingBuffers == 1)
                {
                    // Avoid double loop for the case of a single buffer
                    for (int a = a0; a < a1; a++)
                    {
                        rvec_inc(f[a], fp[0][a]);
                    }
                }
                else
                {
                    for (int a = a0; a < a1; a++)
                    {
                        for (int fb = 0; fb < numContributingBuffers; fb++)
                        {
                            rvec_inc(f[a], fp[fb][a]);
                        }
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

} // namespace

template<typename ForceBufferElementType>
ThreadedForceBuffer<ForceBufferElementType>::ThreadedForceBuffer(const int  numThreads,
                                                                 const bool useEnergyTerms,
                                                                 const int  numEnergyGroups) :
    useEnergyTerms_(useEnergyTerms)
{
    threadForceBuffers_.resize(numThreads);
#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int t = 0; t < numThreads; t++)
    {
        try
        {
            /* Note that thread 0 uses the global fshift and energy arrays,
             * but to keep the code simple, we initialize all data here.
             */
            threadForceBuffers_[t] = std::make_unique<ThreadForceBuffer<ForceBufferElementType>>(
                    t, useEnergyTerms_, numEnergyGroups);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

template<typename ForceBufferElementType>
void ThreadedForceBuffer<ForceBufferElementType>::setupReduction()
{
    const int numBuffers = threadForceBuffers_.size();

    const int numAtoms = threadForceBuffers_[0]->size();

    const int totalNumBlocks =
            (numAtoms + ThreadForceBuffer<ForceBufferElementType>::s_reductionBlockSize - 1)
            >> ThreadForceBuffer<ForceBufferElementType>::s_numReductionBlockBits;

    // Check that all thread buffers have matching sizes
    for (const auto& threadForceBuffer : threadForceBuffers_)
    {
        GMX_RELEASE_ASSERT(threadForceBuffer->size() == numAtoms,
                           "All buffers should have the same size");
        GMX_RELEASE_ASSERT(threadForceBuffer->reductionMask().ssize() == totalNumBlocks,
                           "The block count should match");
    }

    /* Reduce the masks over the threads and determine which blocks
     * we need to reduce over.
     */
    reductionMask_.resize(totalNumBlocks);

    usedBlockIndices_.clear();
    int numBlocksUsed = 0;
    for (int b = 0; b < totalNumBlocks; b++)
    {
        gmx_bitmask_t& mask = reductionMask_[b];

        /* Generate the union over the threads of the bitmask */
        bitmask_clear(&mask);
        for (int t = 0; t < numBuffers; t++)
        {
            bitmask_union(&mask, threadForceBuffers_[t]->reductionMask()[b]);
        }
        if (!bitmask_is_zero(mask))
        {
            usedBlockIndices_.push_back(b);
        }

        if (debug)
        {
            int c = 0;
            for (int t = 0; t < numBuffers; t++)
            {
                if (bitmask_is_set(mask, t))
                {
                    c++;
                }
            }
            numBlocksUsed += c;

            if (gmx_debug_at)
            {
                fprintf(debug, "block %d flags %s count %d\n", b, to_hex_string(mask).c_str(), c);
            }
        }
    }
    if (debug)
    {
        fprintf(debug,
                "Number of %d atom blocks to reduce: %d\n",
                ThreadForceBuffer<ForceBufferElementType>::s_reductionBlockSize,
                int(gmx::ssize(usedBlockIndices_)));
        fprintf(debug,
                "Reduction density %.2f for touched blocks only %.2f\n",
                numBlocksUsed * ThreadForceBuffer<ForceBufferElementType>::s_reductionBlockSize
                        / static_cast<double>(numAtoms),
                numBlocksUsed / static_cast<double>(gmx::ssize(usedBlockIndices_)));
    }
}

template<typename ForceBufferElementType>
void ThreadedForceBuffer<ForceBufferElementType>::reduce(gmx::ForceWithShiftForces* forceWithShiftForces,
                                                         real*                      ener,
                                                         gmx_grppairener_t*         grpp,
                                                         gmx::ArrayRef<real>        dvdl,
                                                         const gmx::StepWorkload& stepWork,
                                                         const int reductionBeginIndex)
{
    if (stepWork.computeForces && !usedBlockIndices_.empty())
    {
        /* Reduce the force buffer */
        GMX_ASSERT(forceWithShiftForces, "Need a valid force buffer for reduction");

        reduceThreadForceBuffers<ForceBufferElementType>(
                forceWithShiftForces->force(), threadForceBuffers_, reductionMask_, usedBlockIndices_);
    }

    const int numBuffers = numThreadBuffers();

    /* When necessary, reduce energy and virial using one thread only */
    if ((stepWork.computeEnergy || stepWork.computeVirial || stepWork.computeDhdl)
        && numBuffers > reductionBeginIndex)
    {
        gmx::ArrayRef<const std::unique_ptr<ThreadForceBuffer<ForceBufferElementType>>> f_t =
                threadForceBuffers_;

        if (stepWork.computeVirial)
        {
            GMX_ASSERT(forceWithShiftForces, "Need a valid force buffer for reduction");

            rvec* gmx_restrict fshift = as_rvec_array(forceWithShiftForces->shiftForces().data());

            for (int i = 0; i < gmx::c_numShiftVectors; i++)
            {
                for (int t = reductionBeginIndex; t < numBuffers; t++)
                {
                    rvec_inc(fshift[i], f_t[t]->shiftForces()[i]);
                }
            }
        }
        if (stepWork.computeEnergy && useEnergyTerms_)
        {
            GMX_ASSERT(ener, "Need a valid energy buffer for reduction");

            for (int i = 0; i < F_NRE; i++)
            {
                for (int t = reductionBeginIndex; t < numBuffers; t++)
                {
                    ener[i] += f_t[t]->energyTerms()[i];
                }
            }
        }

        if (stepWork.computeEnergy)
        {
            GMX_ASSERT(grpp, "Need a valid group pair energy buffer for reduction");

            for (int i = 0; i < static_cast<int>(NonBondedEnergyTerms::Count); i++)
            {
                for (int j = 0; j < f_t[0]->groupPairEnergies().nener; j++)
                {
                    for (int t = reductionBeginIndex; t < numBuffers; t++)
                    {
                        grpp->energyGroupPairTerms[i][j] +=
                                f_t[t]->groupPairEnergies().energyGroupPairTerms[i][j];
                    }
                }
            }
        }
        if (stepWork.computeDhdl)
        {
            GMX_ASSERT(!dvdl.empty(), "Need a valid dV/dl buffer for reduction");

            for (auto i : keysOf(f_t[0]->dvdl()))
            {

                for (int t = reductionBeginIndex; t < numBuffers; t++)
                {
                    dvdl[static_cast<int>(i)] += f_t[t]->dvdl()[i];
                }
            }
        }
    }
}

template class ThreadForceBuffer<RVec>;
template class ThreadedForceBuffer<RVec>;

template class ThreadForceBuffer<rvec4>;
template class ThreadedForceBuffer<rvec4>;

} // namespace gmx
