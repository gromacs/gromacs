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

#ifndef NBLIB_LISTED_FORCES_GPU_PARAM_CUH
#define NBLIB_LISTED_FORCES_GPU_PARAM_CUH

#include <cstdio>

#include <numeric>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "nblib/listed_forces/definitions.h"
#include "nblib/listed_forces/gpu_kernel.cuh"

namespace nblib
{

class ListedInteractionDataGpu
{
public:
    static constexpr int numTypes = TypeListSize<GpuListedTypes>{};

    ListedInteractionDataGpu() = default;

    void update(const ListedInteractionData& cpuParameters)
    {
        auto countInteractions = [this](const auto& interactionElement) {
            using InteractionType = typename std::decay_t<decltype(interactionElement)>::type;
            int typeIndex         = FindIndex<InteractionType, GpuListedTypes>{};

            int numInteractions       = interactionElement.indices.size();
            int numInteractionsPadded = c_threadsPerBlock * iceil(numInteractions, c_threadsPerBlock);

            this->numInteractions_[typeIndex]    = numInteractions;
            this->interactionOffsets_[typeIndex] = numInteractionsPadded;
            this->indexOffsets_[typeIndex] = numInteractionsPadded * (NCenter<InteractionType>{} + 1);
        };

        auto gpuIndices = subsetIndices(GpuListedTypes{}, AllListedTypes{});
        for_each_tuple(countInteractions, tieElements(cpuParameters, gpuIndices));

        numInteractionsTot_ = std::accumulate(interactionOffsets_.begin(), interactionOffsets_.end(), 0);
        std::exclusive_scan(
                interactionOffsets_.begin(), interactionOffsets_.end(), interactionOffsets_.begin(), 0);

        // we store for each thread block which interaction type it will compute
        int                      numBlocks = numInteractionsTot_ / c_threadsPerBlock;
        thrust::host_vector<int> blockTypes(numBlocks);
        for (int typeIndex = 0; typeIndex < numTypes; ++typeIndex)
        {
            int firstBlock = interactionOffsets_[typeIndex] / c_threadsPerBlock;
            int lastBlock  = firstBlock + iceil(numInteractions_[typeIndex], c_threadsPerBlock);
            std::fill(blockTypes.begin() + firstBlock, blockTypes.begin() + lastBlock, typeIndex);
        }
        deviceBlockTypes_ = blockTypes;

        int numIndices = std::accumulate(indexOffsets_.begin(), indexOffsets_.end(), 0);
        std::exclusive_scan(indexOffsets_.begin(), indexOffsets_.end(), indexOffsets_.begin(), 0);

        deviceIndices_.resize(numIndices);

        auto transferType = [this](const auto& interactionElement) {
            using InteractionType   = typename std::decay_t<decltype(interactionElement)>::type;
            constexpr int typeIndex = FindIndex<InteractionType, GpuListedTypes>{};

            // upload parameters to device
            std::get<typeIndex>(this->parametersA_) = interactionElement.parametersA;
            std::get<typeIndex>(this->parametersB_) = interactionElement.parametersB;

            // store pointer in array
            this->deviceParametersA_[typeIndex] =
                    thrust::raw_pointer_cast(std::get<typeIndex>(this->parametersA_).data());
            this->deviceParametersB_[typeIndex] =
                    thrust::raw_pointer_cast(std::get<typeIndex>(this->parametersB_).data());

            int numInteractions = interactionElement.indices.size();
            int typeOffset      = indexOffsets_[typeIndex];
            int indexSize       = NCenter<InteractionType>{} + 1;

            const int* sourceIndices = reinterpret_cast<const int*>(interactionElement.indices.data());
            thrust::copy(sourceIndices,
                         sourceIndices + numInteractions * indexSize,
                         deviceIndices_.begin() + typeOffset);
        };

        for_each_tuple(transferType, tieElements(cpuParameters, gpuIndices));

        devicePointerParametersA_.resize(numTypes);
        devicePointerParametersB_.resize(numTypes);
        devicePointerNumInteractions_.resize(numTypes);
        devicePointerInteractionOffsets_.resize(numTypes);
        devicePointerIndexOffsets_.resize(numTypes);
        thrust::copy(deviceParametersA_.begin(), deviceParametersA_.end(), devicePointerParametersA_.begin());
        thrust::copy(deviceParametersB_.begin(), deviceParametersB_.end(), devicePointerParametersB_.begin());
        thrust::copy(numInteractions_.begin(), numInteractions_.end(), devicePointerNumInteractions_.begin());
        thrust::copy(interactionOffsets_.begin(),
                     interactionOffsets_.end(),
                     devicePointerInteractionOffsets_.begin());
        thrust::copy(indexOffsets_.begin(), indexOffsets_.end(), devicePointerIndexOffsets_.begin());
    }

    const int&  numInteractionsTot() const { return numInteractionsTot_; }
    const auto& numInteractions() const { return numInteractions_; }
    const auto& interactionOffsets() const { return interactionOffsets_; }
    const auto& indexOffsets() const { return indexOffsets_; }

    const auto& parametersA() const { return deviceParametersA_; }
    const auto& parametersB() const { return deviceParametersB_; }

    const int* deviceIndices() const { return thrust::raw_pointer_cast(deviceIndices_.data()); }

    const int* deviceBlockTypes() const
    {
        return thrust::raw_pointer_cast(deviceBlockTypes_.data());
    }

    void* const* deviceParametersA() const
    {
        return thrust::raw_pointer_cast(devicePointerParametersA_.data());
    }
    void* const* deviceParametersB() const
    {
        return thrust::raw_pointer_cast(devicePointerParametersB_.data());
    }
    const int* deviceNumInteractions() const
    {
        return thrust::raw_pointer_cast(devicePointerNumInteractions_.data());
    }
    const int* deviceInteractionOffsets() const
    {
        return thrust::raw_pointer_cast(devicePointerInteractionOffsets_.data());
    }
    const int* deviceIndexOffsets() const
    {
        return thrust::raw_pointer_cast(devicePointerIndexOffsets_.data());
    }

private:
    thrust::device_vector<int> deviceIndices_;

    //! the interaction type of each thread block
    thrust::device_vector<int> deviceBlockTypes_;

    int                        numInteractionsTot_;
    util::array<int, numTypes> numInteractions_;
    //! cumulative sum of number of preceding interactions for each type
    util::array<int, numTypes> interactionOffsets_;

    //! the i-th element of this is the first index of the i-th type in the indices_ array
    util::array<int, numTypes> indexOffsets_;

    util::array<void*, numTypes> deviceParametersA_;
    util::array<void*, numTypes> deviceParametersB_;

    thrust::device_vector<void*> devicePointerParametersA_;
    thrust::device_vector<void*> devicePointerParametersB_;
    thrust::device_vector<int>   devicePointerNumInteractions_;
    thrust::device_vector<int>   devicePointerInteractionOffsets_;
    thrust::device_vector<int>   devicePointerIndexOffsets_;

    using GpuParameterDataList = Map<thrust::device_vector, GpuListedTypes>;
    using GpuParameters        = Reduce<std::tuple, GpuParameterDataList>;

    GpuParameters parametersA_;
    GpuParameters parametersB_;
};

} // namespace nblib

#endif // NBLIB_LISTED_FORCES_GPU_PARAM_CUH
